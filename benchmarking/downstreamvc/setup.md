# Tool Setup

All tools can be installed via conda/mamba into an isolated environment.

## Install conda/mamba

If conda is not installed:

```bash
# Install miniforge (includes mamba)
wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash Miniforge3-Linux-x86_64.sh -b -p $HOME/miniforge3
source $HOME/miniforge3/etc/profile.d/conda.sh
```

## Create the benchmark environment

```bash
mamba create -n varforge-vc -c conda-forge -c bioconda -c defaults \
    bwa-mem2 \
    samtools \
    "gatk4>=4.4" \
    rtg-tools \
    python>=3.10 \
    && conda activate varforge-vc
```

## Prerequisites after environment activation

### Python symlink

The GATK wrapper script calls `python` directly. If only `python3` is available, create a symlink:

```bash
mkdir -p ~/.local/bin && ln -sf $(which python3) ~/.local/bin/python && echo 'export PATH="$HOME/.local/bin:$PATH"' >> ~/.bashrc
source ~/.bashrc
```

### GATK sequence dictionary

Run once before benchmarking. Replace the path with your actual reference location:

```bash
gatk CreateSequenceDictionary -R /path/to/chr22.fa
```

This creates `chr22.dict` alongside the FASTA file. GATK requires it for all pipeline steps.

## Verify installation

```bash
bwa-mem2.avx2 version
samtools --version | head -1
gatk --version
rtg version
```

## Prepare the chr22 reference

Run once before benchmarking:

```bash
cd benchmarking

# BWA-MEM2 index
bwa-mem2.avx2 index chr22.fa

# GATK sequence dictionary
gatk CreateSequenceDictionary -R chr22.fa

# samtools FASTA index (may already exist)
samtools faidx chr22.fa

# RTG-tools SDF (needed for vcfeval)
rtg format -o chr22.sdf chr22.fa
```

## Optional: VarDict

VarDict is a second caller for cross-validation. Install separately as it
requires Java and has a more complex setup:

```bash
mamba install -n varforge-vc -c bioconda vardict-java
```

VarDict is not required to run the benchmark; Mutect2 alone is sufficient.

## Disk space

The full benchmark generates roughly:
- ~2 GB FASTQ per scenario (30x chr22, 150 bp reads)
- ~1.5 GB BAM per scenario (after alignment and sorting)
- ~50 MB VCF files
- Total: ~15 GB for all four scenarios plus intermediate files

Use `--outdir` to direct output to a partition with enough space.
