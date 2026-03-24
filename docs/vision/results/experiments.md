# VarForge: Experiments Before Release

## 1. Mutect2 variant recovery

Simulate known SNVs and indels at four VAF levels: 1%, 5%, 10%, and 50%. Run Mutect2 on the resulting BAM. Measure sensitivity (recall) and precision against the truth VCF using hap.py or vcfeval.

Expected outcome: sensitivity tracks VAF as expected. Precision is high across all levels.

## 2. fgbio UMI pipeline

Simulate duplex UMI data using a known UMI barcode set and family structure. Run fgbio GroupReadsByUmi and CallDuplexConsensusReads. Verify that UMI grouping matches the simulated family structure. Check that duplex consensus calls are correct.

Expected outcome: grouping accuracy matches simulated family sizes. Consensus error rate is within expected range.

## 3. hap.py concordance

Take the truth VCF from experiment 1 and compare it against Mutect2 calls using hap.py. Report SNP and indel sensitivity and precision separately.

## 4. Scaling benchmark

Measure throughput (reads per second) and peak memory at 30x, 100x, and 300x WGS depth on a standard reference (chr1 or whole genome). Run each depth with 1, 4, 8, and 16 threads. Report scaling efficiency.

Expected outcome: memory grows sublinearly with depth. Throughput scales with thread count up to I/O saturation.

## 5. Seed reproducibility

Run the same config with the same seed three times. Verify byte-for-byte identical output. Run with a different seed and verify the output differs.

## 6. cfDNA fragment distribution

Configure a cfDNA run with a known fragment length distribution (e.g., nucleosome-phased 166 bp peak). Measure the actual fragment length distribution in the output BAM using samtools or a custom script. Verify the distribution matches the configured parameters.
