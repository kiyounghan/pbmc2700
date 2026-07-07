#!/bin/bash

# Exit immediately if a command exits with a non-zero status
set -e

echo "=================================================="
echo "   STARTING SINGLE-CELL DATA ACQUISITION PIPELINE "
echo "=================================================="

# --- SECTION 1: QUERY DATASET (2,700 PBMCs from 10x Genomics) ---
echo -e "\n--> Step 1: Verifying Local Query Dataset..."

DATA_URL="https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz"
TARGET_TAR="matrix_2700.tar.gz"
EXTRACT_DIR="filtered_gene_bc_matrices"

# Check if the compressed tarball already exists
if [ -f "$TARGET_TAR" ]; then
    echo "SKIP Tarball '$TARGET_TAR' already exists locally."
else
    echo "Tarball not found. Downloading 2,700 PBMC dataset..."
    curl -o "$TARGET_TAR" "$DATA_URL"
    echo "Download complete."
fi

# Check if the extracted data directory already exists
if [ -d "$EXTRACT_DIR" ]; then
    echo "SKIP Extracted directory '$EXTRACT_DIR/' already exists."
else
    echo "Extracted data directory not found. Unpacking archive..."
    tar -xzf "$TARGET_TAR"
    echo "Extraction complete."
fi


# --- SECTION 2: REFERENCE DATASET (pbmcsca via SeuratData) ---
echo -e "\n--> Step 2: Verifying SeuratData Reference Environment..."

# This block spins up an isolated R instance to safely check the SeuratData cache
R_COMMAND="
if (!requireNamespace('SeuratData', quietly = TRUE)) {
    message('SeuratData package missing. Installing now...')
    install.packages('SeuratData', repos = 'https://cloud.r-project.org')
}

# Load available dataset inventory
installed_data <- SeuratData::AvailableData()

if (!'pbmcsca' %in% installed_data\$Dataset) {
    message('Reference dataset \"pbmcsca\" not found in cache. Downloading via SeuratData...')
    SeuratData::InstallData('pbmcsca')
    message('Reference dataset \"pbmcsca\" successfully downloaded and cached.')
} else {
    message('SKIP Reference dataset \"pbmcsca\" is already cached locally.')
}
"

# Execute the R commands directly from the bash environment
Rscript -e "$R_COMMAND"

echo "======================================================="
echo "  DATASETS Downloaded                                  "
echo "======================================================="
