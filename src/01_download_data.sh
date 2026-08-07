#!/bin/bash

# Exit immediately if a command exits with a non-zero status
set -e

echo "=================================================="
echo "   STARTING SINGLE-CELL DATA ACQUISITION PIPELINE "
echo "=================================================="

# --- SECTION 1: QUERY DATASET (pbmc3k via SeuratData) ---
echo -e "\n--> Step 1: Verifying SeuratData Query Environment..."

# This block spins up an isolated R instance to safely check the SeuratData cache
R_COMMAND="
if (!requireNamespace('SeuratData', quietly = TRUE)) {
    message('SeuratData package missing. Installing now...')
    install.packages('SeuratData', repos = 'https://cloud.r-project.org')
}

# Load available dataset inventory
installed_data <- SeuratData::AvailableData()

if (!'pbmc3k' %in% installed_data\$Dataset) {
    message('Queru dataset \"pbmc3k\" not found in cache. Downloading via SeuratData...')
    SeuratData::InstallData('pbmc3k')
    message('Query dataset \"pbmc3k\" successfully downloaded and cached.')
} else {
    message('SKIP Reference dataset \"pbmc3k\" is already cached locally.')
}
"


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
