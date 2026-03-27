#!/usr/bin/env bash
# Setup script for downstream validation tools (T135, T142).
# Creates a conda environment with all required bioinformatics tools
# and downloads the hg38 chr22 reference.
#
# Usage:
#   bash scripts/setup_downstream.sh
#
# Prerequisites: conda (or mamba) installed.

set -euo pipefail

ENV_NAME="varforge-vc"
REF_DIR="benchmarking/refs"
REF_FA="${REF_DIR}/chr22.fa"

# --------------------------------------------------------------------------
# 1. Create conda environment
# --------------------------------------------------------------------------

if command -v mamba &>/dev/null; then
    SOLVER=mamba
else
    SOLVER=conda
fi

echo "==> Creating conda environment '${ENV_NAME}' with ${SOLVER}..."

${SOLVER} create -y -n "${ENV_NAME}" -c conda-forge -c bioconda -c defaults \
    bwa-mem2 \
    samtools'>=1.17' \
    'gatk4>=4.4' \
    rtg-tools \
    fgbio \
    'python>=3.10'

echo "==> Environment '${ENV_NAME}' created."
echo "    Activate it with:  conda activate ${ENV_NAME}"

# --------------------------------------------------------------------------
# 2. Download hg38 chr22 reference
# --------------------------------------------------------------------------

mkdir -p "${REF_DIR}"

if [ -f "${REF_FA}" ]; then
    echo "==> Reference already exists at ${REF_FA}, skipping download."
else
    echo "==> Downloading hg38 chr22 from UCSC..."
    curl -fSL "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr22.fa.gz" \
        | gunzip > "${REF_FA}"
    echo "==> Downloaded to ${REF_FA}"
fi

# --------------------------------------------------------------------------
# 3. Index the reference (needs the conda env active)
# --------------------------------------------------------------------------

echo ""
echo "==> Next steps (run these after activating the environment):"
echo ""
echo "    conda activate ${ENV_NAME}"
echo ""
echo "    # Index for samtools"
echo "    samtools faidx ${REF_FA}"
echo ""
echo "    # Index for bwa-mem2"
echo "    bwa-mem2.avx2 index ${REF_FA}"
echo ""
echo "    # Sequence dictionary for GATK"
echo "    gatk CreateSequenceDictionary -R ${REF_FA}"
echo ""
echo "    # SDF for rtg vcfeval"
echo "    rtg format -o ${REF_DIR}/chr22.sdf ${REF_FA}"
echo ""
echo "Or run them all at once:"
echo ""
echo "    conda activate ${ENV_NAME} && \\"
echo "    samtools faidx ${REF_FA} && \\"
echo "    bwa-mem2.avx2 index ${REF_FA} && \\"
echo "    gatk CreateSequenceDictionary -R ${REF_FA} && \\"
echo "    rtg format -o ${REF_DIR}/chr22.sdf ${REF_FA}"
echo ""
echo "==> Done. Then run the benchmark with:"
echo "    ./benchmarking/downstreamvc/scripts/run_benchmark.sh \\"
echo "        --ref ${REF_FA} --sdf ${REF_DIR}/chr22.sdf --threads 12"
