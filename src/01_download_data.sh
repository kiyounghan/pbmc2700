#!/bin/bash

# ==============================================================================
# SCRIPT: 01_download_data.sh
# DESCRIPTION: 2 datasets. Downloads and unpacks both the 
#              2,700-cell target query dataset and the 68,000-cell reference atlas.
# ==============================================================================

#  Production error handling
set -e
set -o pipefail

# ==============================================================================
# CONFIGURATION: DIRECTORIIES & URL
# ==============================================================================
# Dataset 1: Target Query (2,700 Cells)
QUERY_DIR="data/pbmc_2700"
QUERY_URL="https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz"
QUERY_TAR="${QUERY_DIR}/pbmc3k_filtered_gene_bc_matrices.tar.gz"

# Dataset 2: Reference Atlas (68,0000 Cells)
REF_DIR="data/reference_atlas"
REF_URL="https://cf.10xgenomics.com/samples/cell-exp/1.1.0/fresh_68k_pbmc_donor_a/fresh_68k_pbmc_donor_a_filtered_gene_bc_matrices.tar.gz"
REF_TAR="${REF_DIR}/fresh_68k_pbmc_donor_a_filtered_gene_bc_matrices.tar.gz"


# ==============================================================================
# PHASE 1: TARGET QUERY (2,700 CELLS)
# ==============================================================================
echo "=== Phase 1: Starting Target Query Pipeline ==="

mkdir -p "$QUERY_DIR"

if [ ! -f "$QUERY_TAR" ]; then
    echo "[INFO] Downloading 2,700 PBMC query data from 10x Genomics..."
    curl -L -o "$QUERY_TAR" "$QUERY_URL"
    echo "[SUCCESS] Query tarball downloaded successfully."
else
    echo "[INFO] Target query tarball already exists. Skipping download."
fi

echo "[INFO] Decompressing query matrix data..."
tar -zxvf "$QUERY_TAR" -C "$QUERY_DIR"

echo "[INFO] Flattening folder structural layout..."
# 10x archives extract into a nested 'filtered_gene_bc_matrices/hg19/' path
if [ -d "${QUERY_DIR}/filtered_gene_bc_matrices/hg19" ]; then
    mv "${QUERY_DIR}/filtered_gene_bc_matrices/hg19/"* "$QUERY_DIR"
    rm -rf "${QUERY_DIR}/filtered_gene_bc_matrices"
fi


# ==============================================================================
#  PHASE 2: REFERENCE ATLAS (68,000 CELLS)
# ==============================================================================
echo "=== Phase 2: Stating Reference Atlas  Pipeline ==="

mkdir -p "$REF_DIR"

if [ ! -f "$REF_TAR" ]; then
    echo "[INFO] Downloading Zheng 68k PBMC reference data from 10x Genomics..."
    curl -L -o "$REF_TAR" "$REF_URL"
    echo "[SUCCESS] Reference atlas tarball downloaded successfully."
else
    echo "[INFO] Reference atlas tarball already exists. Skipping download."
fi

echo "[INFO] Decompressing reference matrix data..."
tar -zxvf "$REF_TAR" -C "$REF_DIR"

echo "[INFO] Flattening folder structural layout..."
# The 68k archive extracts into a nested 'filtered_matrices_mex/hg19/' folder path
if [ -d "${REF_DIR}/filtered_matrices_mex/hg19" ]; then
    mv "${REF_DIR}/filtered_matrices_mex/hg19/"* "$REF_DIR"
    rm -rf "${REF_DIR}/filtered_matrices_mex"
fi


# ==============================================================================
# PIPELINE STATUS CHECK
# ==============================================================================
echo "=== [COMPLETED] Data Staging Complete ==="
echo "[STATUS] Target query dataset ready at: ${QUERY_DIR}/"
echo "[STATUS] Reference atlas dataset ready at: ${REF_DIR}/" 



