#!/bin/bash

# Cleaning
rm -rf data/

# Build directories at the same time
mkdir -p data/pbmc2700

# Download from 10x Genomics
curl -Lo data/pbmc2700/matrix_2700.tar.gz https://cf.10xgenomics.com/samples/cell-exp/1.1.0/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz

# Extract 
tar -xzvf data/pbmc2700/matrix_2700.tar.gz -C data/pbmc2700/
