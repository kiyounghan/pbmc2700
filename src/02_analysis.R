# ==============================================================================
# 1. Libraries & Data
# ==============================================================================
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
if(!require(e1071)) install.packages("e1071", repos="http://cran.us.r-project.org")
library(e1071)

# Ingest sparse count data matrices
data_path <- "data/pbmc2700/filtered_gene_bc_matrices/hg19/"
pbmc.data <- Read10X(data.dir = data_path)
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc2700", min.cells = 3, min.features = 200)

# ==============================================================================
# 2. Data cleaning
# ==============================================================================
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# Global scaling and feature optimization
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
pbmc <- ScaleData(pbmc, features = rownames(pbmc))
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc), verbose = FALSE)

# ==============================================================================
# 3. Clustering
# ==============================================================================
# Extract PCA coordinates
pca_embeddings <- pbmc@reductions$pca@cell.embeddings[, 1:10]

# --- METHOD A: K-Means ---
set.seed(42)
kmeans_output <- kmeans(pca_embeddings, centers = 5)
pbmc$KMeans_Clusters <- as.factor(kmeans_output$cluster)

# --- METHOD B: Hierarchical Clustering  ---
# Expensive method because scales exponentially O(N^2)
distance_matrix <- dist(pca_embeddings)
hierarchical_tree <- hclust(distance_matrix, method = "ward.D2")
pbmc$Hierarchical_Clusters <- as.factor(cutree(hierarchical_tree, k = 5))

# --- METHOD C: Modern Graph-Based Community Detection (Louvain) ---
pbmc <- FindNeighbors(pbmc, dims = 1:10, verbose = FALSE)
pbmc <- FindClusters(pbmc, resolution = 0.5, verbose = FALSE)

# Assign canonical labels to the network graph identities
cell_identities <- c(
  "0" = "CD4+ T Cells", "1" = "CD14+ Monocytes", "2" = "CD4+ T Cells", 
  "3" = "B Cells (MS4A1+)", "4" = "CD8+ T Cells", "5" = "FCGR3A+ Monocytes", 
  "6" = "NK Cells", "7" = "Dendritic Cells"
)
pbmc <- RenameIdents(pbmc, cell_identities)

# ==============================================================================
# 4. Support Vector Machine SVM using KERNEL = Radial Basis Function RBF 
# ==============================================================================
# Mapping cells to using their PCA coordinates
svm_df <- data.frame(PC1 = pca_embeddings[,1], PC2 = pca_embeddings[,2], Label = Idents(pbmc))
svm_model <- svm(Label ~ PC1 + PC2, data = svm_df, kernel = "radial", cost = 1)
pbmc$SVM_Predictions <- predict(svm_model, svm_df)

# ==============================================================================
# 5. Non-Linear Reduction
# ==============================================================================
pbmc <- RunTSNE(pbmc, dims = 1:10)
pbmc <- RunUMAP(pbmc, dims = 1:10, verbose = FALSE)

# ==============================================================================
# 6. VISUALIZATION
# ==============================================================================
# Dashboard 1: Side-by-Side Non-Linear Reduction
tsne_p <- DimPlot(pbmc, reduction = "tsne", label = TRUE, pt.size = 0.4) + ggtitle("t-SNE Projection")
umap_p <- DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.4) + ggtitle("UMAP Projection")
ggsave("results/figures/tsne_vs_umap.png", plot = (tsne_p + umap_p), width = 12, height = 5, dpi = 300)

# Dashboard 2:  Progression of Methods (K-Means vs Hierarchical vs Modern Louvain)
p1 <- DimPlot(pbmc, reduction = "umap", group.by = "KMeans_Clusters", pt.size = 0.3) + ggtitle("Legacy: K-Means")
p2 <- DimPlot(pbmc, reduction = "umap", group.by = "Hierarchical_Clusters", pt.size = 0.3) + ggtitle("Evolutionary: Hierarchical")
p3 <- DimPlot(pbmc, reduction = "umap", pt.size = 0.3) + ggtitle("Modern: Graph (Louvain)")
ggsave("results/figures/clustering_evolution_progression.png", plot = (p1 + p2 + p3), width = 16, height = 5, dpi = 300)

# Dashboard 3: Linear PCA Space vs Supervised SVM Hyperplanes
plot_louvain <- DimPlot(pbmc, reduction = "pca", label = TRUE) + ggtitle("Louvain Graph Classes")
plot_svm     <- DimPlot(pbmc, reduction = "pca", group.by = "SVM_Predictions", label = TRUE) + ggtitle("SVM Decision Hyperplanes (RBF Kernel)")
ggsave("results/figures/svm_boundary_comparison.png", plot = (plot_louvain + plot_svm), width = 12, height = 5, dpi = 300)

# Dashboard 4: Lineage Verification Plots (Canonical Marker Expression)
marker_plot <- FeaturePlot(pbmc, features = c("MS4A1", "CD14", "CD8A", "GNLY"), ncol = 2)
ggsave("results/figures/canonical_markers_verification.png", plot = marker_plot, width = 10, height = 8, dpi = 300)

print("Pipeline operations completed successfully. Matrix figures stored in results/figures/")
