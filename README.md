# For C-SEQTEC Committee's Consideration : 
An end-to-end, reproducible pipeline analyzing 2,700 human Peripheral Blood Mononuclear Cells (PBMCs). This project explores (I) the evolution of unsupervised clustering across simple strategies and non-linear dimensionality reductions, (II) supervised classification using Seurat, (III) a naive Random Forest (RF) model to classify cell lineages across high-dimensional space, and (IV) pseudobulk analysis treating donors as biological conditions.


## Pipeline Architecture
* `src/01_download_data.sh`: Bash script for data download and extraction.
* `src/02_analysis.R`: R script for data transformations and models.
docs/analysis_report.md


## (I) Unsupervised Clustering on pbmc2700 data

## Data Preparation
Data is first filtered using `subset()` then, `NormalizeData()` is applied to correct for differences in sequencing depth (technical variations in total mRNA reads or UMI counts captured per individual cell) so that gene expression levels can be fairly compared across cells. 

`FindVariableFeatures()` is applied to find the top 2,000 highly variable genes (features), hvgs, using Variance-Stabilizing Transformation (vst). These hvgs are converted to a matrix of standardized normal $z$-scores $\mathbf{X} \in \mathbb{R}^{G \times N}$ using `ScaleData()`. 

Finally, `RunPCA()` reads the $X$ and computes Singular Value Decomposition (SVD) $X = U \Sigma V^T$:
* Gene Loadings ($U$) are the "weights" or contribution each gene makes to each Principal Component.
* Singular Values ($\Sigma$) are the values proportional to the variance explained by each Principal Component.
* Cell Embeddings ($V \Sigma$) are the low-dimensional coordinates of each cell in the new low-dimensional PCA space (PC_1, ... , PC_50).   




### 1. The Evolution of Single-Cell Clustering
Compared three distinct clustering methods to evaluate cell lineages:

**K-Means Clustering:** Restricts data points into rigid, spherical clusters around an optimized k=5 centroids by iteratively minimizing the within-cluster sum of squares (WCSS). While computationally simple, it fails to capture non-spherical biological distributions.

* K-Means Top 2 Marker Genes per Classified CellType ranked by `avg_log2FC`
 
![top_KMeans_markers](pbmc2700/results/top_KMeans_markers.png)

**Hierarchical Clustering:** Builds an easy to interpret dendrogram & clustering plot based on pairwise distance matrices (ward.D2: Ward's Minimum Variance Method - True Criterion). While it provides an intuitive visual representation of biological relationships, it scaled poorly with large single-cell datasets (computational complexity: \(\mathcal{O}(N^2)\) and forces cells into strict, irreversible branch splits. Note: k=5 was used to compare against K-Means. 

* Hierarchical Top 2 Marker Genes per Classified CellType ranked by `avg_log2FC`
  
![top_Hierarchical_markers](pbmc2700/results/top_Hierarchical_markers.png) 

**Inter-Clustering Disagreement:** A simple Diagonal Alignment Score = 0.925 confirms strong agreement in grouping between K-Means and Hierarchical clustering is observed. However, a deeper investigation reveals slight mismatch in labeling. After aligning cluster labels via the Hungarian algorithm, Cohen's Unweighted Kappa reached $\kappa$ = 0.889. Meanwhile, Adjusted Rand Index is slightly lower with ARI = 0.819 probably due to how the two algorithms handle the trickier, closely related "T cells" and "NK cells". As seen in the marker genes table, Hierarchical Clustering isolated two canonical hallmark marker genes, `LEF1` and `MAL`, for Cluster 1 - "Naive T Cells" rather than grouping them with other "T cell" subsets. K-Means clustering grouped all "T Cells" (i.e., `TRAT1` and `AQP3`) into a broad Cluster 5 - "T Cells" without separating "Naive T cells". 

The following Cluster Alignment heatmap shows where the labeling discrepancy occurred.

![Cluster_Alignment](pbmc2700/results/Cluster_Alignment.png)

**Inter-Clustering Agreement:** According to Landis & Koch benchmarks, a Kappa value $\kappa$ > 0.8 is considered almost perfect agreement. This almost perfect agreement can be seen from the side-by-side plot below. 

![kmeans + hier](pbmc2700/results/kmeans_hier.png)

**Louvain Clustering:** Builds a Nearest Neighbor network to cluster cells based on graph-based algorithms. This approach scales much more efficiently $(\mathcal{O}(N \cdot k))$ It captures complex non-linear topologies and isolates irregular cell populations without imposing geometric shape constraints. 


![tsne_umap](pbmc2700/results/tsne_umap.png)


### 2. Verification by Canonical Markers

Validate the unsupervised clustering methods and results by generating expression maps for canonical lineage markers. Notably, **`MS4A1`** cleanly marks the **B Cells** cluster, **`CD8A`** marks the **T Cells** cluster since **CD8+ T Cells** is a subpopulation of **T Cells**, **`CD14`** marks Monocytes and Macrophages so it is labeled as **CD14+ Monocytes**, and **`GNLY`** marks **NK Cells** Cluster, confirming the biological accuracy of the groupings.

![canonical_plot](pbmc2700/results/canonical_plot.png)



## (II) Supervised Classification using Seurat {Reference: pbmc16 dataset (used/labeled here as PBMCSCA) & Query: pbmc2700 dataset (used/labeled here as PBMC_2700)}

### 1. Anchor Score Distribution
The PBMCSCA and PBMC_2700 datasets were filtered, then `NormalizeData()`, `FindVariableFeatures()`, `ScaleData()`, and `RunPCA()` were applied. Seurat's supervised label transfer method `FindTransferAnchors()` was applied between the reference and query datasets to find anchors (Mutual Nearest Neighbor cell pairs) to project the reference labels "CellType" directly onto PBMC_2700 cells for classification. There are 4,070 anchors (values between 0 and 1). The summary statistics show median (0.593) and mean (0.596). They are nearly identical. The 3rd Quantile and above are strong alignments with values of 0.778 and above, meaning the distribution is right tail heavy. To be more specific, there is a peak/mode at 1.0 with 459 counts or 11.2% of the total anchors. On the other hand, approximately 3.1% of anchors are valued less than 0.1. These low-valued anchors, between 0 and 0.1, represent ambiguous or noisy cell pairs (e.g., rare or transitional cell states) and soft-weighting is used so that low-scoring anchors contribute very little to final cell-type predictions. Collectively, the distribution of the anchors allow for the expectation of ideal alignment between reference and query because there is high overlap of cell types and low technical noise. The following histogram + summary statistics depicts the distribution of the anchors. 
 
![anchor_histogram](pbmc2700/results/anchor_histogram.png)



### 2. Prediction Score Max
Query label predictions were made using Seurat's `TransferData()`. The prediction score max summary statistics states the 1st Quantile is 0.9723 so this signifies that 75% of cells have a confidence score of 97.23% or higher. In addition, Median = 1.0 and 3rd Quartile = 1.0 so this means at least half of all cells in the query dataset were assigned to a reference cell type with 100% confidence. Across the entire dataset, the average assignment confidence is over 95% because the Mean = 0.9521. Lastly, even the hardest to classify cell in the entire dataset still had a decent ~ 45% maximum probability (Min = 0.4490). Given the prediction score max summary statistics, classification of query cells can be readily accepted. Rather than plotting the distribution or histogram of prediction score max, the following plot visualizes the locations of Prediction Score Max and the summary statistics.   

![pred_score_max_summary_statistics](pbmc2700/results/pred_score_max_summary_statistics.png) 



### 3. Cells of Uncertainty
There are 2,638 cells that have been labeled. A total of 248 cells (~9.40%) fall into the moderate prediction score range (0.50–0.80), while 5 cells (~0.19%) fall into the low prediction score range (0.45–0.50, minimum = 0.449). Of the 253 cells with low-to-moderate probability, 96 (~37.9%) are Cytotoxic T cells (CD8+ T cells).

Cytotoxic T cells simultaneously share a core T-cell lineage with CD4+ T cells, with both expressing identical baseline T-cell markers including CD3D, CD3E, and CD3G. When `FindTransferAnchors()` projects query cells into reference PCA space, these shared T-cell signals can dominate other markers. As a result, cells sitting at the boundary between these clusters split their prediction probabilities across candidate reference labels— allocating weights among Cytotoxic T cells, CD4+ T cells, and NK cells. This probability dispersing directly lowers the maximum prediction score (prediction.score.max).  

Inspecting the prediction score distribution confirms this cross-lineage probability splitting. At high confidence (score $> 0.90$), all 102 cells are cleanly classified as Cytotoxic T cells, whereas 2,263 cells with scores $< 0.10$ are labeled as non-cytotoxic. However, within the intermediate transition zone ($0.10 < \text{score} < 0.90$, $n = 273$), cell assignments split across 129 CD4+ T cells (or their naive, memory, and activated subsets), 121 Cytotoxic T cells, 22 NK cells, and 1 B cell. Notably, CD4+ T-cell subsets form a continuous functional spectrum rather than distinct clusters. This labeling ambiguity between Cytotoxic T cells and NK cells likely stems from their shared cytotoxic effector machinery, as both lineages express identical marker genes such as NKG7, GZMB (Granzyme B), PRF1 (Perforin), and GNLY (Granulysin).

![mid_cells_by_Cyto](pbmc2700/results/mid_cells_by_Cyto.png)


### 4. Predicted CellType by Seurat & Cytotoxic T cells Prediction Score Max Location 
The cells are predicted to be one of nine types. They are labeled and plotted. 

![seurat_predicted_labels](pbmc2700/results/seurat_predicted_labels.png)








## (III) Supervised Classification using Seurat and naive Random Forest {Reference: pbmc16 dataset (used/labeled here as ref) & Query: pbmc2700 dataset (used/labeled here as query)}

### 1. Data Inspection, Analysis, Preparation  

**Data Inspection** Prior to downstream analysis, data inspection and quality control are essential. Inspecting ref@meta.data[1:3, ] reveals several notable columns. For example, the `CellType` column contains 10 unique factor levels, including an "Unassigned" class. The `Cluster` column ranges from 0 to 13. While a `percent.mito` column exists in the reference metadata, it is stored as a fraction; to convert this into a standard percentage scale (0–100%), mitochondrial metrics were recalculated using `ref[["percent.mt"]] <- PercentageFeatureSet(ref, pattern = "^MT-")`. A similar approach was applied to quantify ribosomal read percentages using `ref[["percent.rb"]] <- PercentageFeatureSet(ref, pattern = "^RP[SL]")`.

**Violin Plot, (group.by = "Method"):** Distribution of `nFeature_RNA`, `nCount_RNA`, `percent.mt`, and `percent.rb` stratified by sequencing method. (Cell location was omitted for visual clarity.) Similar plots can be generated by setting `group.by = "CellType"` or `group.by = "orig.ident"`.

The plot illustrates fundamental differences between droplet-based single-cell technologies (e.g., 10x Chromium, Drop-seq, and inDrops) and plate-based, full-length technologies (e.g., Smart-seq2). Droplet methods sequence only the 3' or 5' end of RNA molecules, typically yielding 10k–50k UMI counts per cell. In contrast, full-length methods sequence entire transcripts at much higher coverage per cell, thereby generating 300k–600k+ counts per cell.

![VlnPlot_3_Method](pbmc2700/results/VlnPlot_3_Method.png)

**Multiple Donors:** Multiple Donors: As noted previously, metadata can be grouped by group.by = "orig.ident". The dataset comprises two biological donors: "pbmc1" and "pbmc2". Integrating multiple donors introduces key analytical considerations:
* **Impact on Random Forest Performance:** Incorporating multiple donors allows the Random Forest (RF) model to learn robust cell-type signatures that generalize across individuals by prioritizing universal lineage markers (e.g., CD3D for T cells, MS4A1 for B cells) over donor-specific idiosyncrasies.
* **Risk of Technical Batch Effects:** If "pbmc1" and "pbmc2" were processed on different days or sequenced at varying depths, the classifier risks interpreting technical noise or batch effects as true biological differences.
* **Cell-Type Imbalance Across Donors:** Evaluating lineage distributions between donors revealed significant population shifts. A Chi-Square Test of Independence confirmed that overall cell-type proportions differ significantly between "pbmc1" and "pbmc2" ($\chi^2 = 1774.1$, $\text{df} = 9$, $p < 2.2 \times 10^{-16}$).

To pinpoint which specific lineages drove this imbalance, a pairwise Fisher's Exact Test was performed for each cell type. All cell types exhibited statistically significant donor-level representation differences except for **Plasmacytoid dendritic** cells. To account for this population imbalance and prevent donor-specific overfitting, a stratified 2-fold Leave-One-Group-Out (`orig.ident`) Cross-Validation RF framework was implemented.

The figure below displays a side-by-side bar plot comparing cell-type distributions across donors alongside the corresponding Fisher's Exact Test results:
![Donor_CellType_Distribution](pbmc2700/results/Donor_CellType_Distribution.png)


**Variable genes and Intersection genes:** The top five variable features in the reference dataset are `TAOK1`, `RP11-167N5.5`, `PPBP`, `IGLC2`, and `IGLC3`, whereas the top variable features in the query dataset are `PPBP`, `LYZ`, `S100A9`, `IGLL5`, and `GNLY`.

To train and evaluate the model on a unified feature space, it is critical to intersect feature sets across datasets. Filtering out genes with zero variance yielded **358 common highly variable genes (HVGs)** across 2,282 reference cells and 10 cell types. The top common HVGs include `PPBP`, `PF4`, `PTGDS`, `S100A9`, `S100A8`, and `TUBB1`.


### 2. Modeling: 2-Fold Leave-One-Group-Out (`orig.ident`) Cross-Validation 

The reference dataset exhibits extreme class imbalance, containing 6,872 Cytotoxic T cells compared to only 285 Dendritic cells and 142 Megakaryocytes. Without mitigation, a Random Forest classifier trained on these unweighted distributions will maximize overall accuracy by preferentially predicting majority classes. To ensure rare cell types exert equal influence when determining decision boundaries during tree splits, each cell type was capped at a maximum threshold of 300 cells via stratified subsampling.

**Core Caret Workflow Functions:**
* **`createDataPartition()`:** Performs stratified random sampling, allocating 80% of cells per CellType to the training set.
* **`groupKFold()`:** Generates group-aware fold indices partitioned by donor identity (pbmc1 vs. pbmc2).
* **`trainControl()`:** Encapsulates execution parameters and cross-validation splitting rules into a unified control object.
* **`train()`:** Trains the classifier, evaluates performance across hyperparameter grids, selects the optimal mtry based on Cohen’s Kappa score, and retains the final tuned model.

**Hyperparameter Optimization:**
* **Optimal Parameter:** An optimal `mtry = 24` was selected, achieving a peak Cohen's Kappa score of $\kappa = 0.8472$.
* **Grid Search:** The tuning range (`candidate_mtry <- seq(15, 30, by = 1)`) was refined following multiple coarse parameter sweeps.
![rf_final_plot](pbmc2700/results/rf_final_plot.png)


### 3. Query Projection, Prediction, and Analysis of Fit
Using `FindTransferAnchors()`, 4,070 transfer anchors were identified between the reference and query datasets. The query dataset was subsequently projected into the reference PCA and UMAP embeddings using `RunUMAP()` and `MapQuery()`, and query cell labels were inferred using `predict()`.

**Label Standardization:**
Prior to evaluating classification concordance, string formatting mismatches between predicted labels and original reference metadata were reconciled. Syntactic character substitutions made during Random Forest model training (e.g., replacing spaces and special characters) were standardized back to Seurat's original cell-type naming convention.

Due to these labeling discrepancies, unaligned string names cause `DimPlot()` to map equivalent biological identities to separate color keys, artificially obscuring visual concordance across embeddings.

![p1_p2](pbmc2700/results/p1_p2.png)

Below is a table of the first 9 cells PREDICTED in Seurat's labels, RF labels and RF's cleaned labels. Notice double periods (..) converts to plus signs (+) and single periods (.) changes to spaces ( ) in column 2 and 3, respectively to reconcile string formatting differences and align with Seurat's original metadata conventions.

![Seurat_RF_RF_Clean](pbmc2700/results/Seurat_RF_label.png)

Plot of Seurat's labels vs. RF's cleaned labels. After the name/labeling correction is made, the plot's color scheme is corrected.
![predicted.celltype + rf_model_Labels_clean](pbmc2700/results/predicted_celltype_rf_model_Labels_clean.png) 


Seurat and RF Concordance Heatmap. True Direct Concordance Rate is 92.87% where Concordance Rate is defined by the `mean(Seurat's predicted.celltype == rf_model_Labels_clean * 100)`. 

![Seurat_RF_Concordance](pbmc2700/results/Seurat_RF_Concordance.png)




### 4. Discrepancy analysis

Define a function that yields a Agree/ Disagree values based on  Seurat's `predicted.celltype = rf_model_Labels_clean`. Next, use `FindMarkers()` to find the differentially expressed genes between the Agree and Disagree group of cells. 

Important filter criterias are `p_val_adj` and `avg_log2FC`.

Below is table of **top 5 UP regulated genes** in discrepancy cells where filters : Avg Log2 FC > 0 and filter(p_val_adj < 0.05). 

![top5_UP_filter](pbmc2700/results/top5_UP_filter.png) 
NK cells and cytotoxic/exhausted CD8+ T cells share a substantial portion of their underlying transcriptional program. The primary discrepancy markers identified here— `LAG3`, `GZMK`, `CCL5`, `LYAR`, and `HOPX` —represent key components of this shared effector and exhaustion pathway. This biological overlap explains why automated single-cell classifiers frequently struggle to cleanly resolve these two lineages without targeted sub-clustering or gene-selection filters.

**Calculation of Differential Expression Metrics**

Using the top marker `LAG3` as an exemplar:
* `avg_log2FC`: Log fold change in gene expression, indicating that average expression in the Disagree group is approximately 6.68-fold higher than in the Agree group ($2^{2.74} \approx 6.68$).
* `pct.1` & `pct.2`: The proportion of detected cells (expression > 0) within the Disagree group (`pct.1`) and Agree group (`pct.2`), respectively:
$$\text{pct.1} = \frac{\text{Number of Disagree cells with Gene Expression } > 0}{\text{Total Number of Disagree cells}} = 0.154$$
* `p_val`: Unadjusted p-value is calculated by a non-parametric two-sided Wilcoxon Rank-Sum Test. (At the current moment, trying to understand the statistics employed in `FindMarkers("test.use = "MAST")` Model-based Analysis of Single-cell Transcriptomics. This seems to be more accurate but I am NOT adept in MAST. )
* `p_val_adj`: Adjusted p-value, found by Bonferroni correction, is used to compensate for false positives in multiple testing. It is calculated using p_val_adj = min(1, p_val x N). The formula explains why p_val_adj = 1, in the next table, and therefore, applying filter(p_val_adj < 0.05) is critical.

Below is UP regulated genes with **NO FILTER** applied. Notice how the top gene produces Avg Log2FC = 2.913 but the Adj p-value = 1. Presentation of the following table may seem superfluous and unnecessary but application of filter(p_val_adj < 0.05) is an absolute primary necessity and the potential detrimental impact misleading statistics can have on inferences derived from filter's absence must be highlighted, from the perspective of this statistician.   

![top5_UP_no_filter](pbmc2700/results/top5_UP_no_filter.png)

* Top Fold-Change DOWN regulated Genes with filter and no filter combined.
![Top Fold-Change DOWN regulated Genes](pbmc2700/results/Top_Fold_Change_DOWN_regulated_Genes.png)

### Discrepancy Visualization
* **Spatial Evidence**: The `FeaturePlot()` provided below provides spatial proof on the UMAP that these discrepancy markers are not randomly scattered, but are tightly co-localized to a single, distinct sub-cluster. High expression levels across all four markers concentrate heavily in the lower-middle region of the primary top-right UMAP cluster (chicken drumstick) ($(\text{umap\_1} \in [3, 8], \text{umap\_2} \in [-2, 5])$. This spatial alignment confirms that the shared expression of cytotoxic and exhaustion genes — and the resulting classification disagreement between algorithms — is isolated to a discrete biological subpopulation.

![top_4_discrepancy_genes](pbmc2700/results/top_4_discrepancy_genes.png)


* **Quantifiable Evidence:** The `VlnPlot()` below provides quantitative evidence of the expression variance across the top four markers between cells where classifiers `Agree` versus those where predictions `Disagree`. For `GZMK`, `CCL5`, and `LYAR`, the inner white box plots demonstrate substantial upward shifts in both median expression and interquartile range (IQR) within in the `Disagree` group (spanning normalized expression levels of 2 to 4+), whereas the `Agree` group remains predominantly at zero or low baseline levels. Although, `LAG3` exhibits lower overall baseline expression, it maintains an enriched positive expression (tail density) in the `Disagree` subset. This confirms that the models are not failing due to random technical noise, but rather because these cells actively express a shared cytotoxic, effector, and exhaustion program.

![vlnplot_agree_disagree](pbmc2700/results/vlnplot_agree_disagree.png)


* **Expression Heatmap:** The heatmap shows z-score scaled expression of top marker genes across annotated cell populations. Rows represent individual marker genes and columns represent annotated cell types organized by hierarchical clustering (indicated by the top dendrogram). Color intensity indicates relative gene expression (red = high expression, blue = low expression). Clean diagonal alignment confirms distinct cluster identities, while co-expressed gene bands between `Cytotoxic T cells` and `Natural Killer cells` highlight shared effector transcriptional machinery. 

![Expression_Heatmap](pbmc2700/results/Expression_Heatmap.png)



## (IV) Pseudo-bulk analysis  ("pbmcsca") 

The above methods treated each single cell as an independent statistical observation. However, these methods were based on a flawed statistical assumption: because thousands of cells come from only 2 donor organisms, they are not independent and they share sample-level (`orig.ident`) biological variation. Treating cells as independent replicates inflate sample size and statistical power thereby leading to an unacceptable false-positive rate. Pseudobulk analysis resolves this issue by aggregating transcriptomic counts by sample and cell type into a single bulk profile. This transformation enables established bulk RNA-seq frameworks, such as DESeq2, to properly model sample-to-sample variability and evaluate differential expression at the true biological replicate level.


### scRNA-seq Differential Gene Expression Analysis
Using `FindMarkers()` with filter `group.by = CellType` performs a DE expression analysis and found 5,872 differentially expressed genes between two specific cell populations via a two-sided non-parametric Wilcoxon rank-sum test. The test makes no assumptions about the data following any distribution, making it well-suited for scRNA-seq data, which is heavily zero-inflated and non-normally distributed. However, because it treats individual cells as independent observations, this cell-level approach does not account for biological replicate variability across donors.

For large sample sizes, the distribution of the Mann-Whitney $U$ statistic is asymptotically normal: The distribution of the test statistic $U \sim \mathcal{N}(\mu_U, \sigma_U^2)$ 

Seurat converts $U$ into a standardized $Z$-score: $Z = \frac{U - \mu_U}{\sigma_U}$   


### Pseudo-bulk analysis 

### Data Preparation
* **GetAssayData() & AggregateExpression()**
Raw count data obtained from `GetAssayData(pbmcsca, assay = "RNA", layer = "counts")` are aggregated into pseudobulk profiles using `AggregateExpression(pbmcsca, group.by = c("CellType", "orig.ident"))`. The resulting sparse matrix `dgCMatrix` is generated in Seurat by summing raw transcript counts for each gene across all cells belonging to a given cell-type and sample pair: 
$$\text{Pseudobulk Count}_{\text{GENE, B cell\_pbmc1}} = \sum_{c \in \text{B cells from pbmc1}} \text{Count}_{\text{GENE, cell } c}$$ 

* **Meta-data: df_info**: A metadata data frame (`df_info`) is created from the pseudobulk count matrix so that `DESeq2` can map each aggregated profile to its corresponding biological condition, cell type, and donor (`orig.ident`).

* **Sample Alignment Verification:**
Ensuring that `identical(colnames(pseudobulk_matrix), rownames(df_info))` evaluates to `TRUE` is essential, as DESeq2 requires column names of the count matrix to match the row names of the metadata data frame in both identity and order (`ncol(pseudobulk_matrix) == nrow(df_info)`).

* **Subsetting and Grouping:**
The metadata data frame (`df_info`) is filtered to retain only the cell types and donors of interest. For instance, when comparing `CD4+ T cells` and `Cytotoxic T cells`, the resulting intersect yields four sample profiles: `"CD4+ T cell_pbmc1"`, `"CD4+ T cell_pbmc2"`, `"Cytotoxic T cell_pbmc1"`, and `"Cytotoxic T cell_pbmc2"`.

### DESeqDataSetFromMatrix()
Model Specification and Low-Count Filtering (`DESeqDataSetFromMatrix`)

A `DESeqDataSet` object is constructed from the aggregated pseudobulk count matrix, sample metadata (`df_info`), and experimental design using `DESeqDataSetFromMatrix()`.

* **Design Formula Specification:** The variable order within `design = ~ orig.ident + CellType` is critical. In `DESeq2` design formulas, confounding covariates (such as donor batch effects) are specified first, while the primary variable of interest is placed **last**. This design controls for donor-level variance (`orig.ident`) while evaluating cell-type-specific differential expression (`CellType`).

* **Low-Count Pre-filtering:** Features with fewer than 10 total counts across all combined samples in the subset are filtered out prior to differential testing. Removing low-count genes eliminates non-informative features, reducing the multiple-testing penalty and maximizing statistical power to detect genuinely differentially expressed genes.


### Variance Stabilizing Transformation (`vst()`)
Raw transcript count data are heteroscedastic —where expression variance scales non-linearly as a function of the mean— and thus require transformation prior to downstream visualization and distance calculations. The `vst()` function stabilizes variance across the full dynamic range, yielding homoscedastic values on a $\log_2$-like scale.

Setting `blind = FALSE` ensures that the transformation accounts for the specified experimental design (`design = ~ orig.ident + CellType`). This prevents expected biological variance between cell types and donor batches from being treated as unexplained noise during dispersion trend fitting.     

![var_stab_transf](pbmc2700/results/var_stab_transf.png)

`FGR` expression is roughly at minimum of $2^{(13.066656 - 10.205740)} \approx$ 7.264764-fold higher and at maximum of $2^{(13.587725 - 8.32980)} \approx$ 38.26416-fold higher in `Cytotoxic T cells` compared to `CD4+ T cells`. 


### PCA Plot to visualize donor vs cell-type separation
The following 4 point pseudobulk PCA provides visual verification of the variance splits between cell types and donor batches:

* **PC 1 (47% of variance):** Both `Cytotoxic T cell_pbmc1` and `Cytotoxic T cell_pbmc2` samples (green and purple) sit on the left ($\approx -26$ to $-10$). Both `CD4+ T cell_pbmc1` and `CD4+ T cell_pbmc2` samples (red and cyan) sit on the right ($\approx +11$ to $+25$). This separation confirms that the primary axis of transcriptional variation is driven by lineage differences between cell types.

* **PC 2 (41% of variance):** Both `pbmc1` samples (red and green) are located at the top ($\approx 10$ to $\approx 24$) and both `pbmc2` samples (cyan and purple) are located at the bottom ($\approx -10$ to $\approx -24$). This orthogonal separation confirms that the secondary axis of variation captures donor-specific batch effects.

![PCA_var_stab_transf](pbmc2700/results/PCA_var_stab_transf.png)



### DESeq2

**3 Core Statistical Calculations in `DESeq2`:** 
* 1. Size Factor Estimation (Normalization)
Size factors are estimated using the median-of-ratios method to correct for differences in sequencing depth and prevent compositional bias. The size factors calculated via sizeFactors(DE_seq_object) span a 2.15-fold range across samples:

in $\left( \text{Fold Range} = \frac{\text{Maximum Size Factor}}{\text{Minimum Size Factor}} = \frac{1.640}{0.763} \approx 2.15 \right)$

This fold range indicates that the sample with the highest coverage (`Cytotoxic T cell_pbmc1`) yielded approximately 2.15 times more usable sequencing depth than the sample with the lowest coverage (`CD4+ T cell_pbmc2`).

* 2. Dispersion Estimation 
Mean-Variance Modeling is used to accurately estimate dispersion because pseudobulk counts are overdispersed (variance exceeds the mean).  

* 3. GLM Fitting and Hypothesis Testing (Wald Test)
Negative Binomial Generalized Linear Model is fitted using design formula, calculates Wald test statistics, and derive p-values. `DESeq2` models expected normalized count $\mu_{ij}$ for gene $i$ in sample $j$ using a log link function with an offset for the sample size factor ($s_j$):

$$\log_2(\mu_{ij}) = \beta_0 + \beta_{\text{donor}} \cdot \text{Donor}_j + \beta_{\text{celltype}} \cdot \text{CellType}_j + \log_2(s_j)$$
Model Parameters: $\beta_0$ represents the baseline intercept, while `orig.ident` ($\text{Donor}_j$) and `CellType` ($\text{CellType}_j$) enter the model as dummy-coded covariates.

Effect Size & Testing: The coefficient $\beta_{\text{celltype}}$ directly quantifies the $\log_2$ fold change between cell types while controlling for donor batch variation ($\beta_{\text{donor}}$). The Wald test evaluates the null hypothesis $H_0: \beta_{\text{celltype}} = 0$ by computing $Z = \frac{\beta_{\text{celltype}}}{\text{SE}(\beta_{\text{celltype}})}$, yielding raw $p$-values that are subsequently adjusted for multiple testing using the Benjamini-Hochberg procedure.

**Actuarial Parallel: Negative Binomial GLM in Claim Frequency vs. Pseudobulk Expression**

In actuarial science, insurance claim frequencies ($Y$) are modeled using a Negative Binomial GLM to account for overdispersion ($\alpha > 0$) relative to a standard Poisson model:

$$\mathbb{E}(Y) = \mu = \exp\left(\beta_0 + \sum_{k} \beta_k X_k\right), \quad \text{Var}(Y) = \mu + \alpha \mu^2$$
Where $X_k$ represent policyholder risk factors (e.g., age, geographic location, vehicle class). The quadratic dispersion term $\alpha \mu^2$ prevents insurers from underestimating tail risk and extreme claim events.

In `DESeq2`, transcriptomic read counts are modeled using the exact same statistical framework:

* **Offset Term:** In actuarial GLMs, differences in policy duration are handled via an exposure offset log(Exposure). In `DESeq2`, differences in sequencing depth are handled via the size factor offset $\log_2(s_j)$.

* **Covariates ($X_k$):** Policyholder risk factors correspond to biological condition (`CellType`) and experimental covariates (`orig.ident` / donor batch).

* **Overdispersion ($\alpha$):** In genomics, $\alpha$ captures biological variability and transcriptional bursting across replicates—preventing false positives that would occur if counts were modeled under a naive Poisson assumption ($\text{Var}(Y) = \mu$).



### Effect Size Shrinkage with `lfcShrink(type = "apeglm")`
To obtain accurate, unbiased effect size estimates for ranking genes and downstream visualization, empirical Bayes shrinkage is applied using `lfcShrink()` with `type = "apeglm"`. This method places a heavy-tailed Cauchy prior centered at zero on the $\log_2$ fold changes, shrinking noisy estimates from low-count or high-dispersion genes while preserving large, true biological effects.

* `baseMean`: Average of normalized counts for a single gene calculated across all samples. 
* `log2FoldChange`: The Bayesian posterior mode (peak) of $\log_2$ Fold Change under the Cauchy prior.
* `lfcSE`: Posterior standard deviation of shrunken `log2FoldChange`.
* `stat`: The Wald test statistic calculated as $W = \frac{\log_2(\text{Fold Change})}{\text{lfcSE}}$
* `pvalue`: Unadjusted two sided p-value calculated under $H_0: log_2\ (\text{Fold Change}) = 0$. 
* `padj`: False Discovery Rate (FDR) adjusted $p$-value computed via the Benjamini-Hochberg (BH) procedure.

The `apeglm_df` table is sorted in ascending order by adjusted $p$-value (padj) to prioritize top statistically significant features.
![apeglm_df](pbmc2700/results/apeglm_df.png)


### Heatmap of Top 20 Differentially Expressed Genes (DEGs)

Hierarchical clustering of the top 20 DEGs demonstrates clear transcriptomic separation between cell lineages:

* **Sample Clustering:** The dendrogram resolves into two primary branches, cleanly separating `Cytotoxic T cells` (left) from `CD4+ T cells` (right). This confirms that cell-type-specific biological variation dominates donor-to-donor batch variation.

* **Cytotoxic T-Cell Module:** The lower-left cluster displays marked relative upregulation in `Cytotoxic T cells`. Key features driving this module include canonical cytotoxic effector and chemokine genes: `GZMB`, `GZMH`, `KLRF1`, `CX3CR1`, and `CCL4`.

* **CD4+ T-Cell Module:** The upper-right cluster displays strong relative upregulation in `CD4+ T cells`, highlighted by lineage-associated transcription factors and cell-surface markers: `LEF1`, `MAL`, and `ANKRD55`.

![top20_genes](pbmc2700/results/top20_genes.png)


```text
##################################################################################
NOTE : Volcano map and GO Pathway Enrichment on top DEGs in progress ....
##################################################################################
```


## Local Replication Guidelines

```bash
git clone [https://github.com/kiyounghan/pbmc2700.git](https://github.com/kiyounghan/pbmc2700.git)
cd pbmc2700
bash src/01_download_data.sh
Rscript src/02_analysis.R





