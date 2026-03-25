//! Criterion micro-benchmarks for VarForge hot paths.
//!
//! Run with:
//!   cargo bench
//!
//! Or for a specific benchmark group:
//!   cargo bench -- fragment_sampling

use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};
use rand::rngs::StdRng;
use rand::SeedableRng;
use std::hint::black_box;
use std::io::Write;
use tempfile::TempDir;

use varforge::core::fragment::{CfdnaFragmentSampler, FragmentSampler, NormalFragmentSampler};
use varforge::core::quality::{ParametricQualityModel, QualityModel};
use varforge::core::types::{MutationType, Read, ReadPair};
use varforge::io::fastq::FastqWriter;
use varforge::umi::barcode::{generate_duplex_umi_pair, generate_umi};
use varforge::umi::families::{generate_pcr_copies, UmiFamily};
use varforge::variants::spike_in::{spike_indel, spike_snv};

// ---------------------------------------------------------------------------
// Fragment sampling
// ---------------------------------------------------------------------------

fn bench_fragment_sampling(c: &mut Criterion) {
    let mut group = c.benchmark_group("fragment_sampling");
    group.throughput(Throughput::Elements(1));

    // Normal (WGS) sampler
    let normal_sampler = NormalFragmentSampler::new(300.0, 50.0).unwrap();
    group.bench_function("normal_wgs", |b| {
        let mut rng = StdRng::seed_from_u64(42);
        b.iter(|| black_box(normal_sampler.sample(&mut rng)));
    });

    // cfDNA sampler (more expensive: has ctDNA branch + periodicity)
    let cfdna_sampler = CfdnaFragmentSampler::new(167.0, 334.0, 0.85, 0.05, 20.0, 30.0).unwrap();
    group.bench_function("cfdna", |b| {
        let mut rng = StdRng::seed_from_u64(42);
        b.iter(|| black_box(cfdna_sampler.sample(&mut rng)));
    });

    // Bulk throughput: how many fragments/sec
    for n in [1_000u64, 10_000u64, 100_000u64] {
        group.throughput(Throughput::Elements(n));
        group.bench_with_input(BenchmarkId::new("normal_bulk", n), &n, |b, &n| {
            let mut rng = StdRng::seed_from_u64(42);
            b.iter(|| {
                for _ in 0..n {
                    black_box(normal_sampler.sample(&mut rng));
                }
            });
        });
    }

    group.finish();
}

// ---------------------------------------------------------------------------
// Quality score generation
// ---------------------------------------------------------------------------

fn bench_quality_generation(c: &mut Criterion) {
    let mut group = c.benchmark_group("quality_generation");

    let model = ParametricQualityModel::new(36, 0.003);

    // Per-read cost at common read lengths
    for read_len in [100usize, 150, 250] {
        group.throughput(Throughput::Elements(1));
        group.bench_with_input(
            BenchmarkId::new("generate_qualities", read_len),
            &read_len,
            |b, &len| {
                let mut rng = StdRng::seed_from_u64(42);
                b.iter(|| black_box(model.generate_qualities(len, &mut rng)));
            },
        );
    }

    // Bulk: reads/sec at 150 bp
    for n in [1_000u64, 10_000u64] {
        group.throughput(Throughput::Elements(n));
        group.bench_with_input(BenchmarkId::new("bulk_reads", n), &n, |b, &n| {
            let mut rng = StdRng::seed_from_u64(42);
            b.iter(|| {
                for _ in 0..n {
                    black_box(model.generate_qualities(150, &mut rng));
                }
            });
        });
    }

    group.finish();
}

// ---------------------------------------------------------------------------
// Variant spike-in
// ---------------------------------------------------------------------------

fn bench_variant_spike_in(c: &mut Criterion) {
    let mut group = c.benchmark_group("variant_spike_in");

    // 150 bp read: 37 repeats of "ACGT" (148 bp) + "AC"
    let read_seq: Vec<u8> = b"ACGT".iter().cycle().take(150).copied().collect();
    assert_eq!(read_seq.len(), 150);
    let read_qual = vec![30u8; 150];
    let read_start: u64 = 1000;

    // SNV spike-in — hit case
    let snv_hit = MutationType::Snv {
        pos: 1050,
        ref_base: b'A',
        alt_base: b'T',
    };
    group.throughput(Throughput::Elements(1));
    group.bench_function("snv_hit", |b| {
        b.iter(|| {
            let mut read = Read::new(read_seq.clone(), read_qual.clone());
            black_box(spike_snv(&mut read, read_start, &snv_hit))
        });
    });

    // SNV spike-in — miss case (no overlap)
    let snv_miss = MutationType::Snv {
        pos: 5000,
        ref_base: b'A',
        alt_base: b'T',
    };
    group.bench_function("snv_miss", |b| {
        b.iter(|| {
            let mut read = Read::new(read_seq.clone(), read_qual.clone());
            black_box(spike_snv(&mut read, read_start, &snv_miss))
        });
    });

    // Indel insertion
    let indel_ins = MutationType::Indel {
        pos: 1020,
        ref_seq: vec![b'A'],
        alt_seq: vec![b'A', b'T', b'G'],
    };
    let ref_after = vec![b'N'; 10];
    group.bench_function("indel_insertion", |b| {
        b.iter(|| {
            let mut read = Read::new(read_seq.clone(), read_qual.clone());
            black_box(spike_indel(&mut read, read_start, &indel_ins, &ref_after))
        });
    });

    // Bulk: variants/sec (applying 10 SNVs per read across 1000 reads)
    let snvs: Vec<MutationType> = (0..10)
        .map(|i| MutationType::Snv {
            pos: 1000 + i * 10,
            ref_base: b'A',
            alt_base: b'T',
        })
        .collect();

    group.throughput(Throughput::Elements(1000));
    group.bench_function("bulk_1000_reads_10_snvs", |b| {
        b.iter(|| {
            for _ in 0..1000u32 {
                let mut read = Read::new(read_seq.clone(), read_qual.clone());
                for snv in &snvs {
                    spike_snv(&mut read, read_start, snv);
                }
                black_box(&read);
            }
        });
    });

    group.finish();
}

// ---------------------------------------------------------------------------
// FASTQ writing
// ---------------------------------------------------------------------------

fn make_read_pair(read_len: usize) -> ReadPair {
    ReadPair {
        name: "bench_read".to_string(),
        read1: Read::new(vec![b'A'; read_len], vec![30u8; read_len]),
        read2: Read::new(vec![b'T'; read_len], vec![30u8; read_len]),
        fragment_start: 1000,
        fragment_length: read_len * 2,
        chrom: "chr1".to_string(),
        variant_tags: Vec::new(),
        ref_seq_r1: Vec::new(),
        ref_seq_r2: Vec::new(),
    }
}

fn bench_fastq_writing(c: &mut Criterion) {
    let mut group = c.benchmark_group("fastq_writing");

    // Bytes written per pair: 2 reads × (header ~25 + seq + sep + qual + newlines)
    // At 150 bp: ~2 × (25 + 150 + 2 + 150 + 4) = ~662 bytes/pair, gzip ~2:1 → ~330 bytes
    // We report in terms of pairs written.
    for n_pairs in [100u64, 1_000u64, 10_000u64] {
        let pairs: Vec<ReadPair> = (0..n_pairs as usize).map(|_| make_read_pair(150)).collect();

        group.throughput(Throughput::Elements(n_pairs));
        group.bench_with_input(
            BenchmarkId::new("write_pairs_150bp", n_pairs),
            &pairs,
            |b, pairs| {
                b.iter(|| {
                    let dir = TempDir::new().expect("tempdir");
                    let mut writer = FastqWriter::new(dir.path(), "bench").expect("writer");
                    writer.write_pairs(pairs, "bench").expect("write");
                    writer.finish().expect("finish");
                });
            },
        );
    }

    // Also benchmark writing to an in-memory sink to isolate the format overhead
    // from filesystem + compression overhead.
    group.throughput(Throughput::Elements(1000));
    group.bench_function("write_1000_pairs_in_memory", |b| {
        let pairs: Vec<ReadPair> = (0..1000).map(|_| make_read_pair(150)).collect();
        b.iter(|| {
            let mut buf = Vec::with_capacity(1024 * 1024);
            for (i, pair) in pairs.iter().enumerate() {
                // Write raw FASTQ (uncompressed) to a Vec<u8>
                writeln!(buf, "@bench:{}/1", i + 1).unwrap();
                buf.write_all(&pair.read1.seq).unwrap();
                buf.write_all(b"\n+\n").unwrap();
                let enc1: Vec<u8> = pair.read1.qual.iter().map(|&q| q + 33).collect();
                buf.write_all(&enc1).unwrap();
                buf.write_all(b"\n").unwrap();
                writeln!(buf, "@bench:{}/2", i + 1).unwrap();
                buf.write_all(&pair.read2.seq).unwrap();
                buf.write_all(b"\n+\n").unwrap();
                let enc2: Vec<u8> = pair.read2.qual.iter().map(|&q| q + 33).collect();
                buf.write_all(&enc2).unwrap();
                buf.write_all(b"\n").unwrap();
            }
            black_box(buf)
        });
    });

    group.finish();
}

// ---------------------------------------------------------------------------
// UMI generation and PCR family expansion
// ---------------------------------------------------------------------------

fn bench_umi_generation(c: &mut Criterion) {
    let mut group = c.benchmark_group("umi_generation");

    // Single UMI generation at common lengths
    for umi_len in [8usize, 12, 16] {
        group.throughput(Throughput::Elements(1));
        group.bench_with_input(
            BenchmarkId::new("generate_umi", umi_len),
            &umi_len,
            |b, &len| {
                let mut rng = StdRng::seed_from_u64(42);
                b.iter(|| black_box(generate_umi(len, &mut rng)));
            },
        );
    }

    // Duplex UMI pair
    group.throughput(Throughput::Elements(1));
    group.bench_function("generate_duplex_umi_pair_8bp", |b| {
        let mut rng = StdRng::seed_from_u64(42);
        b.iter(|| black_box(generate_duplex_umi_pair(8, &mut rng)));
    });

    // PCR family expansion
    let original_pair = make_read_pair(150);
    for family_size in [3usize, 10, 50] {
        let family = UmiFamily {
            umi: generate_umi(8, &mut StdRng::seed_from_u64(0)),
            original: original_pair.clone(),
            family_size,
        };
        group.throughput(Throughput::Elements(family_size as u64));
        group.bench_with_input(
            BenchmarkId::new("pcr_family_expand", family_size),
            &family,
            |b, fam| {
                b.iter(|| {
                    let mut rng = StdRng::seed_from_u64(42);
                    black_box(generate_pcr_copies(fam, 0.001, 10, &mut rng))
                });
            },
        );
    }

    // Bulk UMI generation throughput
    for n in [10_000u64, 100_000u64] {
        group.throughput(Throughput::Elements(n));
        group.bench_with_input(BenchmarkId::new("bulk_umi_generation", n), &n, |b, &n| {
            let mut rng = StdRng::seed_from_u64(42);
            b.iter(|| {
                for _ in 0..n {
                    black_box(generate_umi(8, &mut rng));
                }
            });
        });
    }

    group.finish();
}

// ---------------------------------------------------------------------------
// End-to-end simulation pipeline (mini)
// ---------------------------------------------------------------------------
//
// Simulates the inner loop of the simulator: for each fragment, sample size,
// generate reads, apply quality, spike in variants, prepare for writing.

fn bench_pipeline_mini(c: &mut Criterion) {
    let mut group = c.benchmark_group("pipeline_mini");

    let frag_sampler = NormalFragmentSampler::new(300.0, 50.0).unwrap();
    let qual_model = ParametricQualityModel::new(36, 0.003);
    let snv = MutationType::Snv {
        pos: 1050,
        ref_base: b'A',
        alt_base: b'T',
    };

    for n_reads in [1_000u64, 10_000u64] {
        group.throughput(Throughput::Elements(n_reads));
        group.bench_with_input(
            BenchmarkId::new("reads_with_qual_and_snv", n_reads),
            &n_reads,
            |b, &n| {
                b.iter(|| {
                    let mut rng = StdRng::seed_from_u64(42);
                    let mut results = Vec::with_capacity(n as usize);
                    for i in 0..n {
                        let frag_len = frag_sampler.sample(&mut rng);
                        let read_len = (frag_len / 2).min(150);

                        // Generate reference-derived sequence (just A's for bench)
                        let q1 = qual_model.generate_qualities(read_len, &mut rng);
                        let q2 = qual_model.generate_qualities(read_len, &mut rng);

                        // Spike SNV with 50% probability (simulates VAF=0.5)
                        let read_start = 1000u64 + i * 200;
                        let mut r1 = Read::new(vec![b'A'; read_len], q1);
                        spike_snv(&mut r1, read_start, &snv);
                        let r2 = Read::new(vec![b'T'; read_len], q2);

                        results.push(ReadPair {
                            name: format!("read:{}", i),
                            read1: r1,
                            read2: r2,
                            fragment_start: read_start,
                            fragment_length: frag_len,
                            chrom: "chr1".to_string(),
                            variant_tags: Vec::new(),
                            ref_seq_r1: Vec::new(),
                            ref_seq_r2: Vec::new(),
                        });
                    }
                    black_box(results)
                });
            },
        );
    }

    group.finish();
}

// ---------------------------------------------------------------------------
// Criterion configuration
// ---------------------------------------------------------------------------

criterion_group! {
    name = benches;
    config = Criterion::default()
        .sample_size(50)
        .warm_up_time(std::time::Duration::from_secs(1))
        .measurement_time(std::time::Duration::from_secs(3));
    targets =
        bench_fragment_sampling,
        bench_quality_generation,
        bench_variant_spike_in,
        bench_fastq_writing,
        bench_umi_generation,
        bench_pipeline_mini,
}

criterion_main!(benches);
