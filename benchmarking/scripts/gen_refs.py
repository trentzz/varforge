#!/usr/bin/env python3
"""
Generate synthetic reference FASTA files for VarForge benchmarking.

Creates reference genomes of various sizes with realistic ~40% GC content
plus matching .fai index files.
"""

import random
import struct
import sys
import os

def generate_sequence(length: int, gc_fraction: float = 0.40, seed: int = 42) -> bytes:
    """Generate a random DNA sequence with the specified GC content."""
    rng = random.Random(seed)
    at_bases = b'AT'
    gc_bases = b'GC'
    seq = bytearray(length)
    for i in range(length):
        if rng.random() < gc_fraction:
            seq[i] = gc_bases[rng.randint(0, 1)]
        else:
            seq[i] = at_bases[rng.randint(0, 1)]
    return bytes(seq)


def write_fasta_and_fai(path: str, chromosomes: list[tuple[str, int]], gc: float = 0.40):
    """
    Write a FASTA file with the given chromosomes and create a .fai index.

    chromosomes: list of (name, length) tuples
    """
    line_width = 60  # bases per line in the FASTA

    fai_records = []
    offset = 0  # byte offset in the FASTA file

    with open(path, 'wb') as fasta:
        for chrom_idx, (name, length) in enumerate(chromosomes):
            # Write header
            header = f">{name}\n".encode()
            fasta.write(header)
            header_len = len(header)

            # Record where sequence starts
            seq_offset = offset + header_len

            # Generate sequence
            seq = generate_sequence(length, gc, seed=chrom_idx * 1_000_003 + 42)

            # Write sequence in lines of line_width
            full_lines = length // line_width
            remainder = length % line_width

            for i in range(full_lines):
                fasta.write(seq[i * line_width:(i + 1) * line_width])
                fasta.write(b'\n')

            if remainder:
                fasta.write(seq[full_lines * line_width:])
                fasta.write(b'\n')

            # FAI record:
            # NAME  LENGTH  OFFSET  LINEBASES  LINEWIDTH
            # LINEWIDTH includes the newline character
            fai_records.append((name, length, seq_offset, line_width, line_width + 1))

            # Update offset: header + full lines (each line_width + 1 for \n) + remainder line
            offset = seq_offset + full_lines * (line_width + 1)
            if remainder:
                offset += remainder + 1

    # Write .fai index
    fai_path = path + ".fai"
    with open(fai_path, 'w') as fai:
        for name, length, seq_offset, linebases, linewidth in fai_records:
            fai.write(f"{name}\t{length}\t{seq_offset}\t{linebases}\t{linewidth}\n")

    print(f"Written {path} ({os.path.getsize(path):,} bytes) + {fai_path}")


def main():
    out_dir = os.path.dirname(os.path.abspath(__file__))

    # 1 MB reference: ~6 chromosomes of varying sizes summing to ~1 MB
    ref_1mb = [
        ("chr1", 250_000),
        ("chr2", 200_000),
        ("chr3", 180_000),
        ("chr4", 150_000),
        ("chr5", 130_000),
        ("chr6", 100_000),  # total: 1,010,000 ~ 1 MB
    ]

    # 10 MB reference: 8 chromosomes
    ref_10mb = [
        ("chr1", 2_000_000),
        ("chr2", 1_800_000),
        ("chr3", 1_500_000),
        ("chr4", 1_300_000),
        ("chr5", 1_200_000),
        ("chr6", 1_000_000),
        ("chr7",   700_000),
        ("chr8",   600_000),  # total: 10,100,000 ~ 10 MB
    ]

    # 50 MB reference: 12 chromosomes
    ref_50mb = [
        ("chr1",  6_000_000),
        ("chr2",  5_500_000),
        ("chr3",  5_000_000),
        ("chr4",  4_500_000),
        ("chr5",  4_000_000),
        ("chr6",  3_500_000),
        ("chr7",  3_000_000),
        ("chr8",  2_500_000),
        ("chr9",  2_000_000),
        ("chr10", 5_000_000),
        ("chr11", 5_000_000),
        ("chr12", 4_000_000),  # total: 50,500,000 ~ 50 MB
    ]

    write_fasta_and_fai(os.path.join(out_dir, "ref_1mb.fa"), ref_1mb)
    write_fasta_and_fai(os.path.join(out_dir, "ref_10mb.fa"), ref_10mb)
    write_fasta_and_fai(os.path.join(out_dir, "ref_50mb.fa"), ref_50mb)

    print("All reference files generated successfully.")


if __name__ == "__main__":
    main()
