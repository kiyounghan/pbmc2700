# For C-SEQTEC Committee's Consideration : Single-Cell Pipeline: pbmc2700

Query: An end-to-end reproducible pipeline analyzing 2,700 Human Peripheral Blood Mononuclear Cells (PBMCs). This project explores the evolution of models of single-cell clustering strategies, maps 2 non-linear dimensional reductions, and a Random Forest (RF) model to classify cell lineages across high-dimensional spaces.

## Pipeline Architecture
* `src/01_download_data.sh`: Bash script for data downlod and extraction.
* `src/02_analysis.R`: R script for data transformations, various clustering, and ML non-linear boundary models.

## Progression 

### 1. The Evolution of Single-Cell Clustering
Compared three methods to determine clustering efficiency:
* **K-Means:** Restricts data points into rigid, equal-sized spherical clusters. This approach fails to capture continuous developmental trajectories or irregular biological distributions.
* **Hierarchical:** Accurately captures relationships between cells, however there is a problem with a computational complexity from $O(N^2)$ to $O(N^3)$. 
* **Graph-Based (Louvain Modularity):** Builds a Shared Nearest Neighbor (SNN) graph layout that scales linearly. This model efficiently isolates irregular populations without assuming fixed spherical shapes.

![Clustering Evolution Progression](results/clustering_evolution_progression.png)


### 2. Verification by Canonical Markers
Validate the unsupervised clustering methods and results by generating expression maps for canonical lineage markers. Notably, **`MS4A1`** cleanly marks the B-cell cluster, confirming the biological accuracy of our groupings.

![Canonical Marker Profiles](results/canonical_markers_verification.png)

### 3. Supervised Reference Classification (Random Forest - RF)
Used the benchmark training data benchmark from Systematic comparison of single-cell, Nat. Biotech 2020, which profiled a total of 31,021 human PBMCs using 10x Chromium (v2). Trained a 5-fold CV on 80/20 training/test set then fitted RF model to predict Labels of query dataset.  

![Louvain Random Forest](results/Louvain_Random_Forest.png)

### 4. Discrepancy Analysis 
Louvain is not the absolute biological truth. However, use Louvain clusters as an unbiased, unsupervised baseline of the raw data structure. By mapping the supervised Random Forest Labels, the discrepancy analysis allows us to spot exactly where mathematical clustering and supervised biological memory **disagree**. These are the most interesting or ambiguous cell population.

![Discrepancy Random Forest](results/Discrepancies_Random_Forest.png)

## Local Replication Guidelines
```bash
git clone [https://github.com/kiyounghan/pbmc2700.git](https://github.com/kiyounghan/pbmc2700.git)
cd pbmc2700
bash src/01_download_data.sh
Rscript src/02_analysis.R
