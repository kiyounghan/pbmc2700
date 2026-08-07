---
title: "Analysis of Human PBMCs"
author: "Ki Young Han"
date: "2026-06-22"
output: html_document
---


```{r}
rm(list = ls())
```


print("============================================================================================")
print(" I. UNSUPERVISED CLUSTERING (Data = pbmc3k as pbmc2700) ... STARTING                        ")
print("============================================================================================")








print("=======================================================================================") 
print("  0. Libraries... STARTING                                                           ")  
print("=======================================================================================")
```{r setup, include=FALSE}
pacman::p_load(apeglm, BiocManager, caret, clue, DESeq2, doParallel, dplyr, e1071, EnhancedVolcano, fgsea, foreach, ggplot2, gridExtra, gt, irr, kableExtra, knitr, Matrix, mclust, msigdbr, parallel, patchwork, pheatmap, randomForest, ranger, Seurat, SeuratData, tidyverse)
```
print("=======================================================================================")
print("  0. Libraries... COMPLETED                                                          ")
print("=======================================================================================")






print("========================================================================================")
print(" 1. Upload spare data.... STARTING                                                    ") 
print("========================================================================================")

```{r}
pbmc2700 <- LoadData("pbmc3k")
```


```{r}
pbmc2700 <- GetAssayData(pbmc2700, layer = "counts")
pbmc2700 <- CreateSeuratObject(counts = pbmc2700, min.cells = 3, min.features = 200,project = "pbmc2700")
```


```{r}
#str( pbmc2700 , max.level = 20) 
```


# "@i: 29 ..." , "@x:1..." , "@p: [0, 779,2131,..]" so all elements in @i and @x from index 0 to 778 are column/cell 1
```{r}
pbmc2700$RNA$counts[30,1]  # @i: 29 so 29+1 = 30 
```

print("=======================================================================================")
print(" 1. Upload spare data... COMPLETED                                                   ")
print("=======================================================================================")










print("=======================================================================================")
print("  2. Data cleaning... STARTING                                                         ")
print("=======================================================================================")


```{r}
# percentage of mitochondrial reads
pbmc2700[["percent.mt"]] <- PercentageFeatureSet(pbmc2700, pattern = "^MT-")
head(pbmc2700,3)
```


```{r}
colnames(pbmc2700)[1]
```


```{r}
grep("^MT-" , rownames(pbmc2700)) 

sum(pbmc2700$RNA$counts[grep("^MT-" , rownames(pbmc2700)) ,1]) 
pbmc2700$nCount_RNA[1]

aa = sum(pbmc2700$RNA$counts[grep("^MT-" , rownames(pbmc2700)) ,1])  / pbmc2700$nCount_RNA[1] 


message(paste(
  "% mitochondrial reads for", colnames(pbmc2700)[1] , "is:" , round(aa*100, 4), 
  "%. This is equal to percent.mt:",head(pbmc2700)[1,4] )
  )

```


```{r}
pbmc2700[["percent.rb"]] <- PercentageFeatureSet(pbmc2700, pattern = "^RP[SL]")
head(pbmc2700,3)
```


```{r}
head(grep("^RP[SL]" , rownames(pbmc2700))) 

sum(pbmc2700$RNA$counts[grep("^RP[SL]" , rownames(pbmc2700)) ,1]) 
pbmc2700$nCount_RNA[1] 

bb <- sum(pbmc2700$RNA$counts[grep("^RP[SL]" , rownames(pbmc2700)) ,1])  / pbmc2700$nCount_RNA[1] 

# 1057 / 2419 
message(paste(
  "ribosomal % for", colnames(pbmc2700)[1] , "is:" , round(bb*100, 4), 
  "%. This is equal to percent.rb:",head(pbmc2700)[1,5] )
  )
```


```{r}
feat_scat_1 <- FeatureScatter(pbmc2700, feature1 = "percent.rb", feature2 = "percent.mt")

ggsave("/Users/kiyounghan/Desktop/feat_scat_1.png", plot = feat_scat_1 , width  = 12 , height = 10 , dpi = 300 )

feat_scat_1
```


```{r}
ptrb_vlnplot <- VlnPlot(pbmc2700, features = "percent.rb")

ggsave("/Users/kiyounghan/Desktop/ptrb_vlnplot.png", plot = ptrb_vlnplot , width  = 12 , height = 10 , dpi = 300 )

ptrb_vlnplot
```


```{r}
feat_scat_2 <- FeatureScatter(pbmc2700, feature1 = "percent.rb", feature2 = "nCount_RNA")
ggsave("/Users/kiyounghan/Desktop/feat_scat_2.png", plot = feat_scat_2 , width  = 12 , height = 10 , dpi = 300 )

feat_scat_2
```


```{r}

vlnplot1 <- VlnPlot(pbmc2700, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), ncol = 4)
ggsave("/Users/kiyounghan/Desktop/vlnplot1.png", plot = vlnplot1 , width  = 20 , height = 10 , dpi = 300 )

vlnplot1
```

```{r}
summary1 <- summary(pbmc2700@meta.data[, c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb")])
summary1
```


```{r}
pbmc2700 <- subset(pbmc2700, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)  
summary2 <- summary(pbmc2700@meta.data[, c("nFeature_RNA", "nCount_RNA", "percent.mt")])
summary2
```


```{r}

vlnplot2 <- VlnPlot(pbmc2700, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave("/Users/kiyounghan/Desktop/vlnplot2.png", plot = vlnplot2 , width  = 12 , height = 10 , dpi = 300 )
vlnplot2
```

print("=======================================================================================")
print(" 2. Data cleaning... COMPLETED                                                         ")
print("=======================================================================================")





print("=======================================================================================")
print(" 3. Data Normalization, ScaleData, PCA ... STARTING                                    ")
print("=======================================================================================")

```{r}
pbmc2700 <- NormalizeData(pbmc2700, verbose = FALSE)
```


```{r}

head(pbmc2700$RNA$counts@x,10)  
head(pbmc2700$RNA$data@x,10)    
```


```{r}
pbmc2700 <- FindVariableFeatures(pbmc2700, selection.method = "vst", nfeatures = 2000, verbose = FALSE)

top10 <- head(VariableFeatures(pbmc2700), 10)
plot1 <- VariableFeaturePlot(pbmc2700)
LabelPointsplot <- LabelPoints(plot = plot1, points = top10, repel = TRUE) 
LabelPointsplot
ggsave("/Users/kiyounghan/Desktop/LabelPointsplot.png", plot = LabelPointsplot , width  = 12 , height = 10 , dpi = 300 )
```


```{r}
top10
```

# The inputs X
```{r}
pbmc2700 <- ScaleData(pbmc2700, verbose = FALSE)
pbmc2700$RNA$scale.data[1:2,1:3]
```

# The weights W
```{r}
pbmc2700 <- RunPCA(pbmc2700, features = VariableFeatures(object = pbmc2700), verbose = FALSE)
sum( pbmc2700@reductions$pca@feature.loadings[1:2,1:10] != Loadings(pbmc2700, reduction = "pca")[1:2,1:10] )
dim(Embeddings(pbmc2700, reduction = "pca"))
```


print("=======================================================================================")
print(" 3. Data Normalization, ScaleData, PCA ... COMPLETED                                   ")
print("=======================================================================================")








print("=======================================================================================")
print(" 4. Clustering & Dimension Reduction... STARTING                                       ")
print("=======================================================================================")

# --- METHOD A: K-Means --- 

# The Output S = XW
```{r}
pbmc2700_emb <- Embeddings(pbmc2700, reduction = "pca")[, 1:10]
sum( pbmc2700@reductions$pca@cell.embeddings !=  Embeddings(pbmc2700, reduction = "pca") )
```


```{r}
wcss <- vector()
 for (k in 1:10) {
      set.seed(123)
  km <- kmeans(pbmc2700_emb, centers = k, nstart = 10)
  wcss[k] <- km$tot.withinss
}


m <- (wcss[10] - wcss[1]) / (10 - 1)
b <- wcss[1] - m * 1


vertical_heights <- sapply(1:10, 
                           
    function(k) {
  line_y <- m * k + b
  height <- line_y - wcss[k]
  return(height)
                }
  )


optimal_k <- which.max(vertical_heights[1:10])  


data_wcss <- data.frame(
  k = 1:10,
  wcss = wcss,
  line_y = m * (1:10) + b
)

kmeansplot <- ggplot(data_wcss, aes(x = k)) +
  
  geom_line(aes(y = wcss), color = "blue", size = 1) +
  geom_point(aes(y = wcss), color = "black", size = 4) +
  geom_line(aes(y = line_y), color = "gray", linetype = "dashed", size = 1) +
  
  annotate(
    "segment",
    x = optimal_k, 
    y = wcss[optimal_k], 
    xend = optimal_k, 
    yend = m * optimal_k + b,
    color = "red", 
    size = 1,
    linetype = "solid"
  ) +
  

  geom_point(
    data = data_wcss[data_wcss$k == optimal_k, ], 
    aes(y = wcss), 
    color = "red", 
    size = 10, 
    shape = 1, 
    stroke = 1
  ) +
  
  scale_x_continuous(breaks = 1:10) +
  labs(
    x = "Number of Clusters (k)",
    y = "Within-Cluster Sum of Squares (WCSS)",
    title = "Elbow Method (Vertical Distance to Reference Line)",
      ) +
  theme_classic()

kmeansplot
ggsave("/Users/kiyounghan/Desktop/kmeansplot.png", plot = kmeansplot , width  = 12 , height = 10 , dpi = 300 )
```


# K-means clusters
```{r}
set.seed(123)
final_kmeans <- kmeans(pbmc2700_emb, centers = optimal_k, nstart = 10)
pbmc2700$KMeans_Clusters <- as.factor(final_kmeans$cluster)

Idents(pbmc2700) <- "KMeans_Clusters"
KMeans_markers <- FindAllMarkers(pbmc2700, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

top_KMeans_markers <- KMeans_markers %>% 
  group_by(cluster) %>% 
  slice_max(n = 2, order_by = avg_log2FC)

print(top_KMeans_markers)

write.csv(top_KMeans_markers, "top_KMeans_markers.csv", row.names = FALSE)

```

 
```{r}

KMeans_names <- c(
  "1" = "NK Cells",
  "2" = "B Cells",
  "3" = "Platelets",
  "4" = "CD14+ Monocytes" , 
  "5" = "T Cells"
  )


mapped_names <- KMeans_names[as.character(pbmc2700$KMeans_Clusters)]

names(mapped_names) <- colnames(pbmc2700)
pbmc2700$KMeans_cell_types <- mapped_names

Idents(pbmc2700) <- "KMeans_cell_types"

```


```{r}
table(pbmc2700$KMeans_cell_types )
```


```{r}
pbmc2700 <- RunTSNE(pbmc2700, dims = 1:10, verbose = FALSE)
pbmc2700 <- RunUMAP(pbmc2700, dims = 1:10, verbose = FALSE)
print(Reductions(pbmc2700))
```


```{r}
kmeans_legend_order <- c("NK Cells", "B Cells",  "Platelets", "CD14+ Monocytes", "T Cells" )
pbmc2700$KMeans_cell_types <- factor( pbmc2700$KMeans_cell_types, levels = kmeans_legend_order )

Idents(pbmc2700) <- "KMeans_cell_types"


plot_kmeans <- DimPlot(
  pbmc2700, 
  reduction = "umap",  
  pt.size = 1,
  label = TRUE, 
  label.size = 10,
  label.color = "black",
  repel = TRUE
) + 
  ggtitle("K-Means Clustering (k = 5)") + 
  theme(legend.position = "right",
        plot.title = element_text(
          hjust = 0.5,
          size = 30,
          face = "bold",
          color = "black"
          )
        )

plot_kmeans

ggsave("/Users/kiyounghan/Desktop/plot_kmeans.png", plot = plot_kmeans , width  = 15 , height = 15 , dpi = 300 )
```


# --- METHOD B: Hierarchical Clustering  ---
```{r}
library(clue)
library(irr)
# Expensive method because scales exponentially.

set.seed(123)

distance_matrix <- dist(pbmc2700_emb)
hierarchical_tree <- hclust(distance_matrix, method = "ward.D2")
pbmc2700$Hierarchical_Clusters <- as.factor(cutree(hierarchical_tree, k=optimal_k ))     

Idents(pbmc2700) <- "Hierarchical_Clusters"
Hierarchical_markers <- FindAllMarkers(pbmc2700, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

top_Hierarchical_markers <- Hierarchical_markers %>% 
  group_by(cluster) %>% 
  slice_max(n = 2, order_by = avg_log2FC)

print(top_Hierarchical_markers)

write.csv(top_Hierarchical_markers, "top_Hierarchical_markers.csv", row.names = FALSE)

```

    
```{r}

Hierarchical_names <- c(
  "1" = "Naive T Cells",
  "2" = "B Cells",
  "3" = "CD14+ Monocytes",
  "4" = "NK/Cytotoxic Cells",
  "5" = "Platelets"
)

mapped_Hierarchical_names <- Hierarchical_names[as.character(pbmc2700$Hierarchical_Clusters)]

names(mapped_Hierarchical_names) <- colnames(pbmc2700)

pbmc2700$Hierarchical_cell_types <- mapped_Hierarchical_names

Idents(pbmc2700) <- "Hierarchical_cell_types"

```



```{r}    

Hierarchical_legend_order <- c(
  "NK/Cytotoxic Cells",  # Aligns with K-Means "NK Cells"
  "B Cells",             # Matches perfectly
  "Platelets",           # Matches perfectly
  "CD14+ Monocytes" ,     # Matches perfectly
  "Naive T Cells"       # Aligns with K-Means "T Cells"
  )

pbmc2700$Hierarchical_cell_types <- factor(pbmc2700$Hierarchical_cell_types, levels = Hierarchical_legend_order)

Idents(pbmc2700) <- "Hierarchical_cell_types"

plot_hier <- DimPlot(
  pbmc2700, 
  reduction = "umap",         
  pt.size = 1,
  label = TRUE, 
  label.size = 10,
  label.color = "black",
  repel = TRUE
) + 
  ggtitle("Hierarchical Clustering (k = 5)") + 
  theme(legend.position = "right",
        plot.title = element_text(
          hjust = 0.5,
          size = 30,
          face = "bold",
          color = "black"
          )
        )

plot_hier

ggsave("/Users/kiyounghan/Desktop/plot_hier.png", plot = plot_hier , width  = 15 , height = 15 , dpi = 300 )

```

```{r}
table(pbmc2700$Hierarchical_cell_types )
```


```{r}

png("hierarchical_dendrogram.png", width = 11, height = 7, units = "in", res = 300)

par(mar = c(5, 4, 4, 2) + 0.1)

plot(hierarchical_tree, 
     labels = FALSE, 
     main = "Hierarchical Clustering Dendrogram", 
     xlab = "Cells (PBMC 2700)", 
     ylab = "Height (Ward's Linkage)", 
     sub = paste("Cut at k =", length(Hierarchical_names), "Clusters"))

rect_info <- rect.hclust(hierarchical_tree, k = length(Hierarchical_names), border = "red")

box_centers <- sapply(rect_info, function(x) mean(match(x, hierarchical_tree$order)))

cut_height <- mean(hierarchical_tree$height[length(hierarchical_tree$height) - (length(Hierarchical_names)-2): (length(Hierarchical_names)-1)])
y_pos <- cut_height * 1.10  

text(x = box_centers, 
     y = y_pos, 
     labels = Hierarchical_names[as.character(1:length(Hierarchical_names))], 
     col = "darkred", 
     font = 2, 
     cex = 0.85, 
     srt = 45,       
     adj = c(0, 0)) 

dev.off()

```


```{r}
kmeans_hier <- plot_kmeans + plot_hier
kmeans_hier
ggsave("/Users/kiyounghan/Desktop/kmeans_hier.png", plot = kmeans_hier , width  = 24 , height = 15 , dpi = 300 )
```


```{r}

kmeans_hier_table <- table(
  KMeans = pbmc2700$KMeans_cell_types, 
  Hier = pbmc2700$Hierarchical_cell_types )

write.csv(
  as.data.frame.matrix(kmeans_hier_table), 
  file = "kmeans(row)_hier(column)_table.csv", 
  row.names = TRUE )

print(kmeans_hier_table)


# Find optimal label mapping using HUNGARIAN ALGORITHM
mapping <- solve_LSAP(kmeans_hier_table, maximum = TRUE)


new_hierarchical <- factor(pbmc2700$Hierarchical_cell_types)
levels(new_hierarchical) <- levels(factor(pbmc2700$KMeans_cell_types))[mapping]


aligned_kappa <- kappa2(data.frame(
  KMeans = pbmc2700$KMeans_cell_types,
  Hierarchical_Aligned = new_hierarchical), weight = "unweighted")


aa <- (sum(diag(kmeans_hier_table)) / sum(kmeans_hier_table))
adj_Rand_Indx <- adjustedRandIndex(pbmc2700$KMeans_cell_types, pbmc2700$Hierarchical_cell_types)

message(paste("Diagonal Alignment Score:", round(aa, 3)))
message(paste("Cohen's UNWEIGHTED Kappa:", round(aligned_kappa$value, 3)))
message(paste("Adjusted Rand Index:", round( adj_Rand_Indx , 3) ) )

```


```{r}
kmeans_hier_df <- as.data.frame(kmeans_hier_table)

Cluster_Alignment <- ggplot(kmeans_hier_df, aes(x = Hier, y = KMeans, fill = Freq)) +
  geom_tile() +
  geom_text(aes(label = Freq), color = "black") +
  scale_fill_gradient(low = "white", high = "red") +
  theme_minimal() +
  labs(title = "Cluster Alignment: K-Means vs Hierarchical",
       x = "Hierarchical Cell Types",
       y = "K-Means Cell Types")

Cluster_Alignment

ggsave("/Users/kiyounghan/Desktop/Cluster_Alignment.png", plot = Cluster_Alignment , width  = 15 , height = 15 , dpi = 300 )
```


print("===========================================================================================")
print("     Unconstrained Louvain... STARTING                                                     ")
print("===========================================================================================")

# --- METHOD C: Louvain ---
```{r}
set.seed(123)
pbmc2700 <- FindNeighbors(pbmc2700, dims = 1:10)
pbmc2700 <- FindClusters(pbmc2700, resolution = 0.5, algorithm = 1)
pbmc2700$Louvain_Clusters <- pbmc2700$seurat_clusters
```


```{r}
Idents(pbmc2700) <- "Louvain_Clusters"
Louvain_markers <- FindAllMarkers(pbmc2700, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

top_Louvain_markers <- Louvain_markers %>% 
  group_by(cluster) %>% 
  slice_max(n = 2, order_by = avg_log2FC)

print(top_Louvain_markers)

write.csv(top_KMeans_markers, "top_Louvain_markers.csv", row.names = FALSE)

```


```{r}
Louvain_names <- c(
  "0" = "Naive CD4 T Cells",
  "1" = "CD14+ Monocytes",
  "2" = "Memory CD4 T Cells",
  "3" = "B Cells",
  "4" = "CD8 T Cells",
  "5" = "FCGR3A+ Monocytes",
  "6" = "NK Cells",
  "7" = "Dendritic Cells",
  "8" = "Platelets"
)

mapped_louvain <- Louvain_names[as.character(pbmc2700$Louvain_Clusters)]
names(mapped_louvain) <- colnames(pbmc2700)
pbmc2700$Louvain_cell_types <- mapped_louvain
Idents(pbmc2700) <- "Louvain_cell_types"
```


```{r}
kmeans_vs_louvain <- table(
  KMeans = pbmc2700$KMeans_cell_types, 
  Louvain = pbmc2700$Louvain_cell_types )

print(kmeans_vs_louvain)
write.csv(kmeans_vs_louvain, "kmeans(row)_vs_louvain(column).csv", row.names = TRUE)
```

print("========================================================================================")
print(" 4. Clustering & Dimension Reduction.... COMPLETED                                      ")
print("========================================================================================")






print("========================================================================================")
print(" 5. t - SNE , UMAP Visualization.... STARTING                                           ")
print("========================================================================================")

```{r}
tsne <- DimPlot(pbmc2700, 
              reduction = "tsne", 
              pt.size = 0.3,
              label = TRUE, 
              label.size = 7,
              label.color = "black",
              repel = TRUE
              ) + 
  ggtitle("t - SNE") +
  theme(
    legend.position = "right",
        plot.title = element_text(
          hjust = 0.5,
          size = 30,
          face = "bold",
          color = "black"
          )
        )

tsne
ggsave( "/Users/kiyounghan/Desktop/tsne.png" , plot = tsne , width  = 10 , height = 14 , dpi = 300 )
```


```{r}
UMAP <- DimPlot(pbmc2700, 
              reduction = "umap", 
              pt.size = 0.3,
              label = TRUE, 
              label.size = 7,
              label.color = "black",
              repel = TRUE
              ) + 
  ggtitle("UMAP") +
  theme(
    legend.position = "right",
        plot.title = element_text(
          hjust = 0.5,
          size = 30,
          face = "bold",
          color = "black"
          )
        )

UMAP
ggsave( "/Users/kiyounghan/Desktop/UMAP.png" , plot = UMAP , width  = 10 , height = 14 , dpi = 300 )
```


```{r}
tsne_umap <- tsne + UMAP
tsne_umap
ggsave( "/Users/kiyounghan/Desktop/tsne_umap.png" , plot = tsne_umap , width  = 20 , height = 14 , dpi = 300 )

```

print("=========================================================================")
print(" 5. Visualization COMPLETED                                              ")
print("=========================================================================")








print("====================================================================================")
print("  6.  Lineage Verification Plots (Canonical Marker Expression)... STARTING          ")   
print("====================================================================================") 
```{r}
p <-  FeaturePlot(
  pbmc2700,
  features = c("MS4A1", "CD14", "CD8A", "GNLY"),
  ncol = 2,
  reduction = "umap",
  pt.size = 1
) & theme(
    plot.title = element_text(size = 16, face = "bold"), 
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 11)
  ) 

canonical_plot <- p + plot_annotation( title = "Canonical Marker Gene Expression") & 
    theme( plot.title = element_text(size = 22, face = "bold", hjust = 0.5, color = "black") )
  
canonical_plot

ggsave( "/Users/kiyounghan/Desktop/canonical_plot.png" , plot = last_plot(), width = 20 , height = 20, dpi = 300 )

```
print("====================================================================================")
print(" 6.Lineage Verification Plots(Canonical Marker Expression) ...COMPLETED             ")   
print("====================================================================================")   


print("============================================================================================")
print(" I. UNSUPERVISED CLUSTERING (Data = "pbmc3k" used as pbmc2700) ... COMPLETED                ")
print("============================================================================================")












print("========================================================================================================")
print(" I I. SUPERVISED (Reference: PBMCSCA ("pbmcsca") & Query: PBMC_2700 ("pbmc3k")) using Seurat...STARTING ")
print("========================================================================================================")

```{r}
rm(list = ls())
```


```{r setup, include=FALSE}
pacman::p_load(apeglm, BiocManager, caret, clue, DESeq2, doParallel, dplyr, e1071, EnhancedVolcano, fgsea, foreach, ggplot2, gridExtra, gt, irr, kableExtra, knitr, Matrix, mclust, msigdbr, parallel, patchwork, pheatmap, randomForest, ranger, Seurat, SeuratData, tidyverse )
```


# --- QUERY DATA ---
```{r}

PBMC_2700 <- LoadData("pbmc3k")
PBMC_2700 <- GetAssayData(PBMC_2700, layer = "counts")
PBMC_2700 <- CreateSeuratObject(counts = PBMC_2700, min.cells = 3, min.features = 200,project = "PBMC_2700")
PBMC_2700[["percent.mt"]] <- PercentageFeatureSet(PBMC_2700, pattern = "^MT-")
PBMC_2700[["percent.rb"]] <- PercentageFeatureSet(PBMC_2700, pattern = "^RP[SL]")
PBMC_2700 <- subset(PBMC_2700, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5) 
PBMC_2700 <- NormalizeData(PBMC_2700, normalization.method = "LogNormalize", scale.factor = 10000, 
                                         verbose = FALSE)
PBMC_2700 <- FindVariableFeatures(PBMC_2700, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
PBMC_2700 <- ScaleData(PBMC_2700, verbose = FALSE)

PBMC_2700 <- RunPCA(PBMC_2700, features = VariableFeatures(object = PBMC_2700), verbose = FALSE)
PBMC_2700_emb <- Embeddings(PBMC_2700, reduction = "pca")
PBMC_2700 <- RunTSNE(PBMC_2700, dims = 1:10, verbose = FALSE)
PBMC_2700 <- RunUMAP(PBMC_2700, dims = 1:10, verbose = FALSE)

```


# --- REFERENCE DATA ---
```{r}
PBMCSCA <- LoadData("pbmcsca")
```


```{r}
# METADATA ASSIGNMENT: Auto-detect the cell type column name. Initially couldn't find what the name was and had ERROR until this step.
possible_names <- c("cell_type", "celltype", "CellType", "cell.type", "labels")
found_name <- intersect(possible_names, colnames(PBMCSCA@meta.data))

if (length(found_name) > 0) {
  
  PBMCSCA$CellType <- PBMCSCA@meta.data[[found_name[1]]]
  message("Success!!! Finally !!! Found : '", found_name[1], "'")
} else {
  
  print(colnames(PBMCSCA@meta.data))
  stop("Could not find cell type!")
}
```


```{r}
PBMCSCA[["percent.mt"]] <- PercentageFeatureSet(PBMCSCA, pattern = "^MT-")
PBMCSCA[["percent.rb"]] <- PercentageFeatureSet(PBMCSCA, pattern = "^RP[SL]")
```


```{r}
PBMCSCA <- subset(PBMCSCA, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 10 & nCount_RNA < 10000)
PBMCSCA <- NormalizeData(PBMCSCA, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
PBMCSCA <- FindVariableFeatures(PBMCSCA, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
PBMCSCA <- ScaleData(PBMCSCA)
PBMCSCA <- RunPCA(PBMCSCA)
```


```{r}
anchors <- FindTransferAnchors(
  reference = PBMCSCA, 
  query = PBMC_2700, 
  dims = 1:30, 
  reference.reduction = "pca",
  normalization.method = "LogNormalize"
)

anchor_df <- anchors@anchors
anchor_df <- as.data.frame(anchors@anchors)

ref_cells   <- colnames(PBMCSCA)
query_cells <- colnames(PBMC_2700)

anchor_df$PBMCSCA_cell   <- ref_cells[anchor_df$cell1]
anchor_df$PBMC_2700_cell <- query_cells[anchor_df$cell2]

head(anchor_df, 5)
```


```{r}

anchor_df <- as.data.frame(anchors@anchors)

summary(anchor_df$score)

nrow(anchors@anchors)

df <- data.frame(val = anchor_df$score)

Anchor_Score_Distribution <- ggplot(df, aes(x = val)) +
  geom_histogram(bins = 30, fill = "lightblue", color = "black") +
  theme_minimal() +
  labs(title = "Anchor Score Distribution", x = "Anchor Score", y = "Frequency")

ggsave(filename = "/Users/kiyounghan/Desktop/Anchor_Score_Distribution.png", plot = Anchor_Score_Distribution, width = 8, height = 6, dpi= 300)

Anchor_Score_Distribution

sum(anchor_df$score <0.1)
sum(anchor_df$score == 1.0)
```


```{r}

predictions <- TransferData(
  anchorset = anchors, 
  refdata = PBMCSCA$CellType, 
  dims = 1:30
)

PBMC_2700 <- AddMetaData(PBMC_2700, metadata = predictions)


table(PBMC_2700$predicted.id)
sum(table(PBMC_2700$predicted.id))
```


```{r}
Idents(PBMC_2700) <- "predicted.id"
head(PBMC_2700@meta.data[, c("predicted.id", "prediction.score.max")])
```


```{r}

PBMC_2700@meta.data$predicted.id <- predictions$predicted.id
PBMC_2700@meta.data$prediction.score.max <- predictions$prediction.score.max
summary(PBMC_2700$prediction.score.max)

df <- data.frame(val = PBMC_2700@meta.data$prediction.score.max)


prediction_score_max <- ggplot(df, aes(x = val)) +
  geom_histogram(bins = 30, fill = "lightgreen", color = "black") +
  theme_minimal() +
  labs(title = "Prediction Score Max Distribution", x = "Prediction Score Max", y = "Frequency")


ggsave(filename = "/Users/kiyounghan/Desktop/prediction_score_max.png", plot = prediction_score_max, width = 8, height = 6, dpi= 300)

prediction_score_max

sum(PBMC_2700@meta.data$prediction.score.max < 0.8) / 2638   # 254
sum(PBMC_2700@meta.data$prediction.score.max == 1.0) / 2638  # 454
```


```{r}
dd = (sum(PBMC_2700@meta.data$prediction.score.max > 0.8) / length(PBMC_2700@meta.data$prediction.score.max) )* 100
     
message( round(dd,6), "% of cells are lableled with a prediction score max of 0.8 or greater.")  
```


```{r}
df <- data.frame(score = PBMC_2700$prediction.score.max)

label_transfer_confidence_dist <- ggplot(df, aes(x = score)) +
  geom_histogram(binwidth = 0.05, fill = "lightblue", color = "black", boundary = 0) +
  geom_vline(xintercept = 0.8, color = "red", linetype = "dashed", linewidth = 1) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 2500)) +
  labs(
    title = "Labels Transfer Confidence Distribution (red dashed line = 0.800)",
    subtitle = "90.371494% of cells are labeled with a prediction score max of 0.800 or greater",
    x = "Max Prediction Score",
    y = "Frequency"
  ) +
  theme_classic()

label_transfer_confidence_dist

ggsave("/Users/kiyounghan/Desktop/label_transfer_confidence_dist.png", plot = label_transfer_confidence_dist, width = 8, height = 6, dpi = 300)
```


```{r}
unique(PBMC_2700@meta.data$predicted.id)
```

 
```{r}

seurat_predicted_labels <- DimPlot(PBMC_2700, group.by = "predicted.id", label = TRUE, label.box = TRUE,repel = TRUE, label.size = 4,label.color = "black") +
  ggtitle("Query (PBMC_2700) CellType Predictions")

seurat_predicted_labels

ggsave("/Users/kiyounghan/Desktop/seurat_predicted_labels.png" , plot = seurat_predicted_labels , width  = 20 , height = 14 , dpi = 300 )
```
 

```{r}

mod_cells <- PBMC_2700@meta.data[ 
  PBMC_2700$prediction.score.max > 0.5 & PBMC_2700$prediction.score.max < 0.80 , 
  c("predicted.id", "prediction.score.max") ]

head(mod_cells)
nrow(mod_cells)
length(PBMC_2700$prediction.score.max)
nrow(mod_cells) / length(PBMC_2700$prediction.score.max)
```


```{r}
low_cells <- PBMC_2700@meta.data[ 
  PBMC_2700$prediction.score.max < 0.50 & 
  PBMC_2700$prediction.score.max > min(PBMC_2700$prediction.score.max) , 
  c("predicted.id", "prediction.score.max") ]

head(low_cells)
nrow(low_cells)
length(PBMC_2700$prediction.score.max)
nrow(low_cells) / length(PBMC_2700$prediction.score.max)
```


```{r}
PBMC_2700@meta.data
```


```{r}

table(PBMC_2700$predicted.id[PBMC_2700$prediction.score.max < 0.80])
sum(table(PBMC_2700$predicted.id[PBMC_2700$prediction.score.max < 0.80]))
```


```{r}
library(tibble)
library(readr)

low_mod_idx <- which(PBMC_2700$prediction.score.max  < 0.80 & 
                      PBMC_2700$prediction.score.Cytotoxic.T.cell > 0  & 
                     PBMC_2700$prediction.score.CD4..T.cell  >0 & 
                      PBMC_2700$prediction.score.Natural.killer.cell >0
                     )

NK_CD4_split_CytotoxicT <- head(PBMC_2700@meta.data[low_mod_idx, c("predicted.id" , "prediction.score.Cytotoxic.T.cell" , "prediction.score.CD4..T.cell" , "prediction.score.Natural.killer.cell" ) ] , 8 )
NK_CD4_split_CyT_row <- data.frame( cell_id = rownames(NK_CD4_split_CytotoxicT), NK_CD4_split_CytotoxicT, row.names = NULL )
NK_CD4_split_CyT_row
write_csv(NK_CD4_split_CyT_row , "/Users/kiyounghan/Desktop/NK_CD4_split_CyT_row.csv")
```


```{r}

head(PBMC_2700@meta.data[low_mod_idx, grep("prediction.score", colnames(PBMC_2700@meta.data) ) ] )
```


```{r}
summary(PBMC_2700@meta.data$prediction.score.Cytotoxic.T.cell)
sum(PBMC_2700@meta.data$prediction.score.Cytotoxic.T.cell < 0.1)
length(PBMC_2700@meta.data$prediction.score.Cytotoxic.T.cell)
sum(PBMC_2700@meta.data$prediction.score.Cytotoxic.T.cell < 0.1) / length(PBMC_2700@meta.data$prediction.score.Cytotoxic.T.cell)

sum(PBMC_2700@meta.data$predicted.id == "Cytotoxic T cell" & PBMC_2700@meta.data$prediction.score.Cytotoxic.T.cell < 0.1)   # number of cells labeled as Cytotoxic T cell when prob < 0.1
sum(PBMC_2700@meta.data$predicted.id == "Cytotoxic T cell")
sum(PBMC_2700@meta.data$predicted.id != "Cytotoxic T cell")
```


```{r}
df <- data.frame(score = PBMC_2700@meta.data$prediction.score.Cytotoxic.T.cell)

Cytotoxic_T_Prediction_Score <- ggplot(df, aes(x = score)) +
  geom_histogram(binwidth = 0.05, fill = "lightblue", color = "black", boundary = 0) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 2500)) +
  labs(
    title = "Cytotoxic T cell Prediction Score",
    subtitle = "2263 are labeled as non-Cytotoxic T cells & 223 labeled as Cytotoxic T cells",
    x = "Prediction Score",
    y = "Frequency"
  ) +
  theme_classic()

Cytotoxic_T_Prediction_Score

ggsave("/Users/kiyounghan/Desktop/Cytotoxic_T_Prediction_Score.png", plot = Cytotoxic_T_Prediction_Score, width = 8, height = 6, dpi = 300)
```


```{r}
high_cells <- PBMC_2700@meta.data[PBMC_2700$prediction.score.Cytotoxic.T.cell > 0.9 , ] 
table(high_cells$predicted.id)
prob_CyT_2 = sum(high_cells$predicted.id == "Cytotoxic T cell") / (sum(table(high_cells$predicted.id))) * 100
message("Given Cytotoxic T cell Probability Score > 0.9, ", round(prob_CyT_2,4), " % of cells are labeled as Cytotoxic T cell.")
```


```{r}
mid_cells <- PBMC_2700@meta.data[PBMC_2700$prediction.score.Cytotoxic.T.cell > 0.1 & PBMC_2700$prediction.score.Cytotoxic.T.cell < 0.9 , ] 
table(mid_cells$predicted.id)
#sum(table(mid_cells$predicted.id))
prob_CyT_3 = sum(mid_cells$predicted.id == "Cytotoxic T cell") / (sum(table(mid_cells$predicted.id))) * 100
message("There are a total of ", sum(table(mid_cells$predicted.id)), " cells with 0.1 < Prediction Score Cytotoxic T cell < 0.9.")
message("Given 0.1 < Cytotoxic T cell Probability Score < 0.9, ", round(prob_CyT_3,4), " % of cells are labeled as Cytotoxic T cell.")
```


```{r}

mid_cells <- PBMC_2700@meta.data[
  PBMC_2700$prediction.score.Cytotoxic.T.cell > 0.1 & 
  PBMC_2700$prediction.score.Cytotoxic.T.cell < 0.9, ] 

table_data <- as.data.frame(table(Predicted_ID = mid_cells$predicted.id))

table_data <- table_data[order(-table_data$Freq), ]

table_grob <- tableGrob(
  table_data,
  rows = NULL,
  cols = c("Predicted Cell Type", "Cell Count"),
  theme = ttheme_minimal(
    base_size = 10,
    core = list(bg_params = list(fill = c("grey", "white"))),
    colhead = list(fg_params = list(fontface = "bold"))))

df <- data.frame(score = PBMC_2700@meta.data$prediction.score.Cytotoxic.T.cell)

Cytotoxic_T_Prediction_Score <- ggplot(df, aes(x = score)) +
  geom_histogram(binwidth = 0.05, fill = "lightblue", color = "black", boundary = 0) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 2500)) +
  labs(
    title = "Cytotoxic T Cell Prediction Score",
    subtitle = "2263 labeled as non-Cytotoxic T & 223 labeled as Cytotoxic T",
    x = "Prediction Score",
    y = "Frequency"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 9))

mid_cells_by_Cyto <- Cytotoxic_T_Prediction_Score + wrap_elements(table_grob) + plot_layout(widths = c(2, 1))

mid_cells_by_Cyto

ggsave("/Users/kiyounghan/Desktop/mid_cells_by_Cyto.png", plot = mid_cells_by_Cyto, width = 11, height = 5, dpi = 300)
```


```{r}
low_cells <- PBMC_2700@meta.data[PBMC_2700$prediction.score.Cytotoxic.T.cell <= 0.1 , ] 
table(low_cells$predicted.id)
prob_CyT = ((sum(low_cells$predicted.id == "CD4+ T cell") + sum(low_cells$predicted.id == "Natural killer cell"))/ (sum(table(low_cells$predicted.id)))) * 100
message("Given Cytotoxic T cell Probability Score < 0.1, ", round(prob_CyT,2), " % of cells are labeled as CD4+ T cell or Natural killer cell.")
```


```{r}
low_conf_barcodes <- rownames(PBMC_2700@meta.data[PBMC_2700$prediction.score.max < 0.80, ])

head(low_conf_barcodes)
```


```{r}
PBMC_2700$low_confidence <- PBMC_2700$prediction.score.max < 0.80

DimPlot(PBMC_2700, group.by = "low_confidence", cols = c("grey", "red")) +
  ggtitle("Query Prediction Confidence Score < 0.80 (Red)") +
  scale_color_manual(
    labels = c("High Confidence (>= 0.80)", "Low Confidence (< 0.80)"),
    values = c("grey", "red")
  )
```


```{r}
Prediction_Score_Max_Location <- FeaturePlot(PBMC_2700, features = "prediction.score.max") +
  scale_color_viridis_c(
    option = "inferno",  # Options: "magma", "plasma", "inferno", "viridis", "cividis", "rocket" , "mako"
    direction = 1,
    limits = c(0, 1)
  ) +
  ggtitle("Query (PBMC_2700) Prediction Score Max Location") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

Prediction_Score_Max_Location

ggsave( "/Users/kiyounghan/Desktop/Prediction_Score_Max_Location.png" , plot = Prediction_Score_Max_Location , width  = 20 , height = 14 , dpi = 300 )
```


```{r}
sum_vals <- summary(PBMC_2700$prediction.score.max)
sum_df <- data.frame(
  Metric = names(sum_vals),
  Value  = round(as.numeric(sum_vals), 4)
)

large_table <- ttheme_default(
  core = list(
    fg_params = list(fontsize = 22, fontface = "plain"), 
    bg_params = list(fill = c("grey", "white"))       
  ),
  colhead = list(
    fg_params = list(fontsize = 22, fontface = "bold"), 
    bg_params = list(fill = "grey")
  ),
  padding = unit(c(8, 8), "mm")                           
)

table_plot <- tableGrob(sum_df, rows = NULL, theme = large_table)


Prediction_Score_Max_Location <- FeaturePlot(PBMC_2700, features = "prediction.score.max") +
  scale_color_viridis_c(
    option = "inferno",
    direction = 1,
    limits = c(0, 1)
  ) +
  ggtitle("Query (PBMC_2700) Prediction Score Max Location") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))


pred_score_max_summary_statistics <- wrap_elements(table_plot) + Prediction_Score_Max_Location +
  plot_layout(widths = c(1, 2.5))  


pred_score_max_summary_statistics
ggsave("/Users/kiyounghan/Desktop/pred_score_max_summary_statistics.png" , plot = pred_score_max_summary_statistics , width  = 20 , height = 14 , dpi = 300 )
```


```{r}
seurat_and_CyTo_Prediction <- seurat_predicted_labels + Prediction_Score_Max_Location
seurat_and_CyTo_Prediction
ggsave("/Users/kiyounghan/Desktop/seurat_and_CyTo_Prediction.png" , plot = seurat_and_CyTo_Prediction , width  = 20 , height = 14 , dpi = 300 )

```


```{r}

CyTo_prediction <- mid_cells_by_Cyto + Cytotoxic_T_Prediction_Score + plot_layout(widths = c(1.2, 1))
CyTo_prediction
ggsave("/Users/kiyounghan/Desktop/CyTo_prediction.png" , plot = CyTo_prediction , width  = 20 , height = 14 , dpi = 300 )
```

print("=============================================================================================================")
print(" II.SUPERVISED (Reference: PBMCSCA ("pbmcsca") & Query: PBMC_2700 ("pbmc3k")) using Seurat...COMPLETED       ")
print("=============================================================================================================")







print("==========================================================================================================")
print(" III.SUPERVISED (Reference: ref ("pbmcsca") & Query: query ("pbmc3k")) by naive Random Forest...STARTING  ")
print("==========================================================================================================")

# --- REFERENCE DATA ---
```{r}
ref <- LoadData("pbmcsca")  # DON'T DO SEURAT OBJECTS. ALREADY DONE by package
```


```{r}
possible_names <- c("cell_type", "celltype", "CellType", "cell.type", "labels")
found_name <- intersect(possible_names, colnames(ref@meta.data))

if (length(found_name) > 0) {
  ref$CellType <- ref@meta.data[[found_name[1]]]
  message("Success!!! Found and copied cell labels from metadata column: '", found_name[1], "'")
} else {
  print(colnames(ref@meta.data))
  stop("Could not find cell type metadata. Look at the printed column names above and choose the right one!")
}
```


```{r}
ref@meta.data[1:3, ]
```


```{r}
#str(ref,20)
unique(ref@meta.data$Method)
```


```{r}
unique(ref@meta.data$CellType)
sort(as.numeric(unique(ref@meta.data$Cluster)))
```


```{r}
ref@meta.data[ref$Cluster == "1" , ][37:38,]
```


```{r}
ref[["percent.mt"]] <- PercentageFeatureSet(ref, pattern = "^MT-")
ref[["percent.rb"]] <- PercentageFeatureSet(ref, pattern = "^RP[SL]")
```


```{r}
mt_genes <- grep(pattern = "^MT-", x = rownames(ref), value = TRUE)
mt_genes
```


```{r}
mt_counts <- GetAssayData(ref, layer = "counts")[mt_genes, ]
mt_counts[1:3, 1:5]
```


```{r}

VlnPlot_1_CellType <- VlnPlot(ref, features = c("nFeature_RNA", "nCount_RNA"), pt.size = 0, group.by = "CellType", ncol=2) & 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)  )

VlnPlot_1_CellType
ggsave( "/Users/kiyounghan/Desktop/VlnPlot_1_CellType.png" , plot = VlnPlot_1_CellType , width  = 20 , height = 14 , dpi = 300 )

############################################################################################################################################

VlnPlot_2_CellType <- VlnPlot(ref,features = c("percent.mt","percent.rb"), pt.size = 0, group.by = "CellType", ncol=2) & 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)  )

VlnPlot_2_CellType 
ggsave( "/Users/kiyounghan/Desktop/VlnPlot_2_CellType.png" , plot = VlnPlot_2_CellType , width  = 20 , height = 14 , dpi = 300 )

#############################################################################################################################################
```


```{r}
VlnPlot_3_CellType <- VlnPlot(ref,features = c("nFeature_RNA", "nCount_RNA","percent.mt","percent.rb"), pt.size = 0, group.by = "CellType", ncol=4) & 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)  )

VlnPlot_3_CellType 
ggsave( "/Users/kiyounghan/Desktop/VlnPlot_3_CellType.png" , plot = VlnPlot_3_CellType , width  = 20 , height = 14 , dpi = 300 )
```


```{r}

VlnPlot_1_Method <- VlnPlot(ref,features = c("nFeature_RNA", "nCount_RNA"), pt.size = 0, group.by = "Method",ncol = 2) & 
   theme(axis.text.x = element_text(angle = 45, hjust = 1)  )

VlnPlot_1_Method 
ggsave( "/Users/kiyounghan/Desktop/VlnPlot_1_Method .png" , plot = VlnPlot_1_Method  , width  = 20 , height = 14 , dpi = 300 )

########################################################################################################################

VlnPlot_2_Method <- VlnPlot(ref, features = c("percent.mt","percent.rb"), pt.size = 0, group.by = "Method", ncol = 2) & 
    theme( axis.text.x = element_text(angle = 45, hjust = 1)  )

VlnPlot_2_Method 
ggsave( "/Users/kiyounghan/Desktop/VlnPlot_2_Method .png" , plot = VlnPlot_2_Method  , width  = 20 , height = 14 , dpi = 300 )

########################################################################################################################
```


```{r}
VlnPlot_3_Method <- VlnPlot(ref, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), pt.size = 0, group.by = "Method", ncol = 4) & 
    theme( axis.text.x = element_text(angle = 45, hjust = 1)  )

VlnPlot_3_Method 
ggsave( "/Users/kiyounghan/Desktop/VlnPlot_3_Method .png" , plot = VlnPlot_3_Method  , width  = 20 , height = 14 , dpi = 300 )
```


```{r}
levels(ref$orig.ident)
```
    
    
```{r}

VlnPlot_1_orig_ident <- VlnPlot(ref,features = c("nFeature_RNA", "nCount_RNA"), pt.size = 0, group.by = "orig.ident",ncol = 2) & 
   theme(axis.text.x = element_text(angle = 45, hjust = 1)  )

VlnPlot_1_orig_ident 
ggsave( "/Users/kiyounghan/Desktop/VlnPlot_1_orig_ident .png" , plot = VlnPlot_1_orig_ident  , width  = 20 , height = 14 , dpi = 300 )

########################################################################################################################

VlnPlot_2_orig_ident <- VlnPlot(ref, features = c("percent.mt","percent.rb"), pt.size = 0, group.by = "orig.ident", ncol = 2) & 
    theme( axis.text.x = element_text(angle = 45, hjust = 1)  )

VlnPlot_2_orig_ident 
ggsave( "/Users/kiyounghan/Desktop/VlnPlot_2_orig_ident.png" , plot = VlnPlot_2_orig_ident  , width  = 20 , height = 14 , dpi = 300 )

########################################################################################################################
```


```{r}
VlnPlot_3_orig_ident <- VlnPlot(ref, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rb"), pt.size = 0, group.by = "orig.ident", ncol = 4) & 
    theme( axis.text.x = element_text(angle = 45, hjust = 1)  )

VlnPlot_3_orig_ident 
ggsave( "/Users/kiyounghan/Desktop/VlnPlot_3_orig_ident.png" , plot = VlnPlot_3_orig_ident  , width  = 20 , height = 14 , dpi = 300 )

```   

  
```{r}
cont_table <- table(CellType = ref$CellType, Original_Identity = ref$orig.ident)
prop_table <- prop.table(cont_table, margin = 2) * 100

cont_prop <- as.data.frame(cbind(cont_table,round(prop_table, 2)), )
colnames(cont_prop) <-c("Cell Counts pbmc1",
                        "Cell Counts pbmc2",
                        "Cell Types Percentage in pbmc1" ,
                        "Cell Types Percentage in pbmc2"
                        )

cont_prop

cont_prop %>%
  tibble::rownames_to_column(var = "Cell Type") %>%
  write_csv("/Users/kiyounghan/Desktop/cont_prop.csv")
```


```{r}

df_plot <- as.data.frame(prop_table)
colnames(df_plot) <- c("CellType", "Donor", "Percentage")

Cell_Type_Distribution_Across_Donor <- ggplot(df_plot, aes(x = Donor, y = Percentage, fill = CellType)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_classic() +
  labs(
    title = "Cell Type Distribution Across Donor Batches",
    x = "Donor Batch",
    y = "Percentage of Total Cells (%)"
  ) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

Cell_Type_Distribution_Across_Donor

ggsave( "/Users/kiyounghan/Desktop/Cell Type Distribution Across Donor Batches.png" , plot = Cell_Type_Distribution_Across_Donor , width  = 20 , height = 14 , dpi = 300 )
```   


```{r}

chisq_res <- chisq.test(cont_table)
print(chisq_res)
```


```{r}

col_totals <- colSums(cont_table)

per_cell_results <- do.call(rbind, lapply(rownames(cont_table), function(ct) {
  
  count_ct <- cont_table[ct, ]
  count_others <- colSums(cont_table) - count_ct
  mat_2x2 <- rbind(count_ct, count_others)
  
  test <- fisher.test(mat_2x2)
  
  pct_pbmc1 <- round((count_ct[1] / col_totals[1]) * 100, 2)
  pct_pbmc2 <- round((count_ct[2] / col_totals[2]) * 100, 2)
  
  data.frame(
    CellType = ct,
    pbmc1_Count = count_ct[1],
    pbmc2_Count = count_ct[2],
    pbmc1_Pct = pct_pbmc1,
    pbmc2_Pct = pct_pbmc2,
    p_value = test$p.value,
    stringsAsFactors = FALSE
  )
}))


per_cell_results$adj_p_value <- p.adjust(per_cell_results$p_value, method = "BH")
per_cell_results$Significant <- ifelse(per_cell_results$adj_p_value < 0.05, "Yes", "No")

per_cell_table <- per_cell_results[, c("CellType", "pbmc1_Count", "pbmc2_Count", "pbmc1_Pct", "pbmc2_Pct", "p_value", "adj_p_value", "Significant")]

write.csv(per_cell_table, file = "per_cell_table.csv", row.names = TRUE)

print(per_cell_table)
```



```{r}

display_table <- per_cell_table
display_table$p_value <- formatC(display_table$p_value, format = "e", digits = 2)
display_table$adj_p_value <- formatC(display_table$adj_p_value, format = "e", digits = 2)

stats_grob <- tableGrob(
  display_table,
  rows = NULL,
  cols = c("Cell Type", "PBMC1 N", "PBMC2 N", "PBMC1 %", "PBMC2 %", "p-value", "FDR adj p", "Sig"),
  theme = ttheme_minimal(
    base_size = 8, 
    core = list(
      bg_params = list(fill = c("grey", "white"), col = NA) # col = NA prevents opaque border boxes
    ),
    colhead = list(fg_params = list(fontface = "bold"))  )  )


Cell_Type_Distribution_Across_Donor <- ggplot(df_plot, aes(x = Donor, y = Percentage, fill = CellType)) +
  geom_bar(stat = "identity", position = "stack", width = 0.6) +
  theme_classic() +
  labs(
    title = "Cell Type Distribution Across Donors",
    x = "Donor Batch",
    y = "Percentage of Total Cells (%)"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 11),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 8)
  ) +
  guides(fill = guide_legend(ncol = 2)) 

side_by_side_donor <- Cell_Type_Distribution_Across_Donor + 
  wrap_elements(panel = stats_grob) + 
  plot_layout(widths = c(1, 1.8)) 

side_by_side_donor

ggsave("/Users/kiyounghan/Desktop/Donor_CellType_Distribution.png", plot = side_by_side_donor, width = 15, height = 6, dpi = 300 )
```


```{r}
ref <- subset(ref, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 10 & nCount_RNA < 10000)

ref <- NormalizeData(ref, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
ref <- FindVariableFeatures(ref, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
ref <- ScaleData(ref, verbose = FALSE)
ref <- RunPCA(ref, features = VariableFeatures(ref), verbose = FALSE)
```


```{r}
ref10 <- head(VariableFeatures(ref),10)
ref10
plotref10 <- VariableFeaturePlot(ref)
plotref10a <- LabelPoints(plot = plotref10 , points = ref10 , repel=TRUE )
plotref10a

ggsave( "/Users/kiyounghan/Desktop/plotref10a.png" , plot = plotref10a , width  = 20 , height = 14 , dpi = 300 )
```


```{r}

dim(Embeddings(ref, reduction = "pca"))
Embeddings(ref, reduction = "pca")[1:2,1:5]
```


```{r}
dim(Loadings(ref, reduction = "pca"))
Loadings(ref, reduction = "pca")[1:2,1:5]
```


```{r}
Stdev(ref, reduction = "pca")
```


# --- QUERY DATA ---
```{r}
query <- LoadData("pbmc3k")
query <- GetAssayData(query, layer = "counts")

query <- CreateSeuratObject(counts = query, min.cells = 3, min.features = 200,project = "pbmc3k")

query[["percent.mt"]] <- PercentageFeatureSet(query, pattern = "^MT-")
query[["percent.rb"]] <- PercentageFeatureSet(query, pattern = "^RP[SL]")
query <- subset(query, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5) 

query <- NormalizeData(query, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
query <- FindVariableFeatures(query, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
query <- ScaleData(query, verbose = FALSE)
query <- RunPCA(query, features = VariableFeatures(query), verbose = FALSE)
query <- RunTSNE(query, dims = 1:10, verbose = FALSE)
query <- RunUMAP(query, dims = 1:10, verbose = FALSE)
```


```{r}
ref10 <- head(VariableFeatures(query),10)
ref10
plotquery10 <- VariableFeaturePlot(query)
plotquery10a <- LabelPoints(plot = plotquery10 , points = ref10 , repel=TRUE )
plotquery10a

ggsave("/Users/kiyounghan/Desktop/plotquery10a.png" , plot = plotquery10a , width  = 20 , height = 14 , dpi = 300 )

```


```{r}

dim(Embeddings(query, reduction = "pca"))
Embeddings(query, reduction = "pca")[1:2,1:5]
```


```{r}
dim(Loadings(query, reduction = "pca"))
Loadings(query, reduction = "pca")[1:2,1:5]
```


```{r}
Stdev(query, reduction = "pca")
```


print("====================================================================================")
print(" 5. Intersection genes  ... STARTING                                                ")
print("====================================================================================")
```{r}
ref_genes <- VariableFeatures(ref)
query_genes <- VariableFeatures(query)

length(ref_genes) 
length(query_genes)

head(ref_genes,5)
head(query_genes,5)

hvgs <- intersect(ref_genes, query_genes)

#seems redundant but no error after this line
final_features <- intersect(hvgs, intersect(rownames(ref), rownames(query)) )  

head(final_features)
length(final_features)
```

print("====================================================================================")
print(" 5. Intersection genes...COMPLETED                                                  ")
print("====================================================================================")




print("====================================================================================")
print(" 6. Data management for Seurat, caret, & RF compatibility...STARTING                ")
print("====================================================================================")

```{r}

X_ref   <- t(GetAssayData(ref, layer = "data")[final_features, ])
X_query <- t(GetAssayData(query, layer = "data")[final_features, ])

colnames(X_ref)   <- make.names(colnames(X_ref))
colnames(X_query) <- make.names(colnames(X_query))

gene_vars <- apply(as.matrix(X_ref), 2, var)
varying_genes <- names(gene_vars[gene_vars > 0])

message("Removed ", length(colnames(X_ref)) - length(varying_genes), " zero-variance gene(s).")

X_ref   <- as.matrix(X_ref[, varying_genes])
X_query <- as.matrix(X_query[, varying_genes])

X_ref[is.na(X_ref)] <- 0
X_query[is.na(X_query)] <- 0

stopifnot(all(colnames(X_ref) == colnames(X_query)))

message("Ready for Random Forest: ", ncol(X_ref), " genes across ", nrow(X_ref), " reference cells.")
```

print("====================================================================================")
print(" 6. Data management for Seurat, caret, & RF compatibility...COMPLETED               ")
print("====================================================================================")








print("====================================================================================")
print(" 7.  80 % Train / 20 % Test STRATIFIED!! Sampling...STARTING                        ")
print("====================================================================================")

```{r}
table(ref@meta.data$CellType)
```


```{r}
set.seed(123)

train_idx <- createDataPartition(ref$CellType, p = 0.80, list = FALSE)[, 1] # [,1] extracts it as a clean integer VECTOR


train_labels <- ref$CellType[train_idx]    

downsampled_sub_idx <- unlist(
  lapply(
    split(seq_along(train_labels), train_labels), 
    function(idx) { if (length(idx) > 300) sample(idx, 300) else idx }
          )
                              )

final_train_idx <- train_idx[downsampled_sub_idx]


test_idx <- setdiff(seq_len(nrow(X_ref)), train_idx)

X_train <- as.matrix(X_ref[final_train_idx, ])
X_test<- as.matrix(X_ref[test_idx, ])

all_celltypes <- make.names(unique(ref$CellType))
y_train <- factor(make.names(ref$CellType[final_train_idx]), levels = all_celltypes)
y_test <- factor(make.names(ref$CellType[test_idx]), levels = all_celltypes)

donor_train <- ref$orig.ident[final_train_idx] 

train_df <- as.data.frame(X_train)
train_df$Label <- y_train
train_df$Donor <- donor_train 
```

print("====================================================================================")
print(" 7.  80 % Train / 20 % Test STRATIFIED!! Sampling...COMPLETED                       ")
print("====================================================================================")







print("====================================================================================")
print(" 8. Parallel Processors ...STARTING                                                 ")            
print("====================================================================================")

```{r}

donor_folds <- caret::groupKFold(train_df$Donor, k = length(unique(train_df$Donor)))

cv_control <- caret::trainControl(
  method      = "cv",
  index     = donor_folds,
  savePredictions = "final",
  allowParallel = TRUE,
  classProbs  = TRUE)

max_feats <- ncol(X_train)
candidate_mtry <- seq(15, 30, by = 1)
valid_mtry <- candidate_mtry[candidate_mtry <= max_feats]

rf_grid <- expand.grid(
  mtry          = valid_mtry,
  splitrule     = "gini",
  min.node.size = 1)

cores_use <- 10
cl <- makeCluster(cores_use)
registerDoParallel(cl)

set.seed(123)
rf_model <- caret::train(
  Label ~ .,
  data  = train_df[, !colnames(train_df) %in% "Donor"], 
  method  = "ranger",
  trControl = cv_control,
  tuneGrid  = rf_grid,
  num.trees  = 250,
  metric  = "Kappa",
  importance = "impurity")

test_df <- as.data.frame(X_test)
test_predictions <- predict(rf_model, newdata = test_df)

try(stopCluster(cl), silent = TRUE)
foreach::registerDoSEQ()
rm(cl)
gc()

best_results <- merge(rf_model$bestTune, rf_model$results)
best_results
```

print("====================================================================================")
print(" 8. Parallel Processors ...COMPLETED                                               ")
print("====================================================================================")





print("====================================================================================")
print(" 9. RF Model Plot...STARTING                                                        ")
print("====================================================================================")

```{r}

rf_model_final_plot <- ggplot(rf_model) +
  geom_point(size = 4, color = "blue") +
  labs(
    title = "mtry (at 24) vs Kappa Statistic (max = 0.8472009)",
    x = "mtry #",
    y = "Kappa Statistic"
  ) +
  theme_minimal(base_size = 12) + 
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"))

rf_model_final_plot
ggsave( "/Users/kiyounghan/Desktop/rf_model_final_plot.png" , plot = rf_model_final_plot , width  = 20 , height = 14 , dpi = 300 )
```

print("====================================================================================")
print(" 9. RF Model Plot...COMPLETED                                                        ")
print("====================================================================================")


print("====================================================================================")
print(" 10. RF Model Accuracy on Test Set...STARTING                                       ")
print("====================================================================================")

# NOT USEFUL !!! NEED TO CLEAN UP LABELS 
```{r}
conf_matrix <- confusionMatrix(test_predictions, y_test)

cat("\n====================================================================================\n")
cat(sprintf("True Validation Accuracy: %.2f%%\n", conf_matrix$overall["Accuracy"] * 100))
cat(sprintf("True Kappa Statistic: %.4f\n", conf_matrix$overall["Kappa"]))
cat("====================================================================================\n\n")

by_class_df <- as.data.frame(conf_matrix$byClass[, c("Sensitivity", "Specificity", "Balanced Accuracy")])
write.csv(by_class_df, file = "RF_per_class_metrics.csv", row.names = TRUE)


print(conf_matrix$byClass[, c("Sensitivity", "Specificity", "Balanced Accuracy")])
```

print("====================================================================================")
print(" 9. RF Model Accuracy on Test Set...COMPLETED                                       ")
print("====================================================================================")




print("====================================================================================")
print(" 10. Query Projection & Prediction...STARTING                                       ")
print("====================================================================================")

```{r}
anchors <- FindTransferAnchors(
  reference            = ref,
  query                = query,
  dims                 = 1:30,
  reference.reduction  = "pca",
  normalization.method = "LogNormalize"
)

anchors_df <- anchors@anchors
anchors_df <- as.data.frame(anchors@anchors)

ref_cells   <- colnames(ref)
query_cells <- colnames(query)

anchors_df$ref_cell   <- ref_cells[anchors_df$cell1]
anchors_df$query_cell <- query_cells[anchors_df$cell2]
```


```{r}
head(anchors_df,5)
```


```{r}
stats <- summary(anchors_df$score)
summary_text <- paste(names(stats), format(round(stats, 3), nsmall = 3), sep = ": ", collapse = "\n")

anchor_histogram <- ggplot(anchors_df, aes(x = score)) +
  geom_histogram(binwidth = 0.05, fill = "lightblue", color = "black") +
  theme_classic() +
  labs(title = "4070 Anchors Score Distribution", x = "Anchor Score", y = "Frequency") +
  annotate(
    "label", 
    x = 0.02, y = Inf,            
    label = summary_text, 
    hjust = 0, vjust = 1.1, 
    size = 3.5, 
    family = "mono" )

anchor_histogram

ggsave("/Users/kiyounghan/Desktop/anchor_histogram.png", plot = anchor_histogram, width = 11, height = 5, dpi = 300)
```


```{r}
ref <- RunUMAP(ref, dims = 1:30, return.model = TRUE)

query <- MapQuery(
  anchorset             = anchors,
  reference             = ref,
  query                 = query,
  refdata               = list(celltype = "CellType"), 
  reference.reduction   = "pca",
  reduction.model       = "umap"
)
```


```{r}

query_preds <- predict(rf_model, newdata = as.data.frame(X_query))

query$rf_model_Labels <- factor(query_preds)

```


```{r}

"rf_model_Labels" %in% colnames(query@meta.data)
```


print("==============================================================================================")
print(" 11. Cross Evaluation: Random Forest vs Seurat ... STARTING                                   ")
print("==============================================================================================")

```{r}

raw_seurat <- as.character(query$predicted.celltype)
raw_rf     <- as.character(query$rf_model_Labels)

head(data.frame(Seurat = raw_seurat, RF = raw_rf), 10)
```


# mismatched name - formatting . MOST IMPORTANT SECTION. 
```{r}
# Clean up Random Forest labels to match Seurat's original formatting !!!! This caused all the troubles.
rf_clean <- query$rf_model_Labels %>%
  as.character() %>%
  gsub("\\.\\.", "+ ", .) %>%  # Convert ".." back to "+ "
  gsub("\\.", " ", .)         # Convert remaining "." back to spaces


Seurat_RF_label <- data.frame(Seurat = as.character(query$predicted.celltype), RF = raw_rf, RF_Clean = rf_clean)
head(Seurat_RF_label)

write.csv(Seurat_RF_label, file = "/Users/kiyounghan/Desktop/Seurat_RF_label_predictions.csv", row.names = FALSE)


query$rf_model_Labels_clean <- factor(rf_clean)

seurat_clean <- as.character(query$predicted.celltype)
true_agreement <- mean(seurat_clean == rf_clean, na.rm = TRUE) * 100

cat("\n====================================================================\n")
cat(sprintf("TRUE Direct Concordance Rate: %.2f%%\n", true_agreement))
cat("====================================================================\n\n")

cross_tab_clean <- table(Seurat = seurat_clean, Random_Forest = rf_clean)
print(cross_tab_clean)
```




```{r}

p1 <- DimPlot(query, reduction = "ref.umap", group.by = "predicted.celltype", label = TRUE, repel = TRUE) +
  ggtitle("Query projected on Ref UMAP (Seurat Transferred)") +
  theme_minimal() +
  theme(legend.position = "none", plot.title = element_text(face = "bold", hjust = 0.5))

p2 <- DimPlot(query, reduction = "ref.umap", group.by = "rf_model_Labels", label = TRUE, repel = TRUE) +
  ggtitle("Query projected on Ref UMAP (Random Forest Labels)") +
  theme_minimal() +
  theme(legend.position = "none", plot.title = element_text(face = "bold", hjust = 0.5))

p1 | p2
```


# group.by = "rf_model_Labels"
```{r}
library(patchwork)

p1 <- DimPlot(query, reduction = "umap", group.by = "predicted.celltype", label = TRUE, repel = TRUE) +
  ggtitle("Query projected on UMAP (Seurat Transferred)") +
  theme_minimal() +
  theme(legend.position = "none", plot.title = element_text(face = "bold", hjust = 0.5))

p2 <- DimPlot(query, reduction = "umap", group.by = "rf_model_Labels", label = TRUE, repel = TRUE) +
  ggtitle("Query projected on UMAP (Random Forest Labels)") +
  theme_minimal() +
  theme(legend.position = "none", plot.title = element_text(face = "bold", hjust = 0.5))


p1_p2 <- p1+p2
p1_p2
ggsave( "/Users/kiyounghan/Desktop/p1_p2.png" , plot = p1_p2 , width  = 20 , height = 14 , dpi = 300 )
```


# group.by = "rf_model_Labels_clean"
```{r}

predicted.celltype <- DimPlot(query, reduction = "umap", group.by = "predicted.celltype", label = TRUE, repel = TRUE) +
  ggtitle("Query projected on Ref UMAP (Seurat Transferred)") +
  theme_minimal() +
  theme(legend.position = "none", plot.title = element_text(face = "bold", hjust = 0.5))

rf_model_Labels_clean <- DimPlot(query, reduction = "umap", group.by = "rf_model_Labels_clean", label = TRUE, repel = TRUE) +
  ggtitle("Query projected on Ref UMAP (Random Forest Labels)") +
  theme_minimal() +
  theme(legend.position = "none", plot.title = element_text(face = "bold", hjust = 0.5))

p0 <- predicted.celltype + rf_model_Labels_clean

p0

ggsave( "/Users/kiyounghan/Desktop/predicted.celltype + rf_model_Labels_clean.png", plot = p0, width  = 20 , height = 15 , dpi = 400)
```


```{r}
cross_tab_clean <- table(
  Seurat = query$predicted.celltype, 
  Random_Forest = query$rf_model_Labels_clean)

prop_table <- prop.table(cross_tab_clean, margin = 1)

Seurat_RF_Concordance <- pheatmap(
  prop_table,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = colorRampPalette(c("white", "grey", "lightblue"))(100),
  main = "Seurat vs. Random Forest Concordance (92.87% Overall Agreement)",
  display_numbers = TRUE,
  number_format = "%.2f",
  fontsize_number = 9,
  angle_col    = 45)

Seurat_RF_Concordance

ggsave( "/Users/kiyounghan/Desktop/Seurat_RF_Concordance.png", plot = Seurat_RF_Concordance, width  = 20 , height = 15 , dpi = 400)
```


```{r}
query$model_agreement <- ifelse(
  as.character(query$predicted.celltype) == as.character(query$rf_model_Labels_clean),
  "Agree",
  "Disagree"  )

Cell_Annotation_Diagreements <- DimPlot(
  query, 
  reduction = "umap", 
  group.by = "model_agreement", 
  cols = c("Agree" = "grey", "Disagree" = "red"),
  pt.size = 0.6
) +
  ggtitle("Spatial Distribution of Cell Annotation Diagreements (~10%)") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5) )

Cell_Annotation_Diagreements

ggsave( "/Users/kiyounghan/Desktop/Spatial Distribution of Cell Annotation Diagreements.png", plot = Cell_Annotation_Diagreements, width  = 20 , height = 15 , dpi = 400)
```


```{r}

Idents(query) <- "model_agreement"

discrepancy_up_DEGs <- FindMarkers(
  query,
  ident.1 = "Disagree",
  ident.2 = "Agree",
  only.pos = TRUE,         
  min.pct = 0.10,          
  logfc.threshold = 0.25)

top10_discrepancy_genes <- discrepancy_up_DEGs %>%
  filter(p_val_adj < 0.05) %>%
  arrange(p_val_adj, desc(avg_log2FC)) %>%
  head(10)

cat("\n*** TOP 10 UPREGULATED DISCREPANCY GENES ***\n")
print(top10_discrepancy_genes)

```


```{r}
discrepancy_up_DEGs <- FindMarkers(
  query,
  ident.1 = "Disagree",
  ident.2 = "Agree",
  only.pos = TRUE,         
  min.pct = 0.10,          
  logfc.threshold = 0.25)

top5_UP_filter <- discrepancy_up_DEGs %>%
  tibble::rownames_to_column(var = "Gene") %>%
  filter(p_val_adj < 0.05) %>%
  arrange(p_val_adj, desc(avg_log2FC)) %>%
  head(5) %>%
  select(Gene, avg_log2FC, pct.1, pct.2, p_val,p_val_adj)

top5_UP_filter %>%
  kable(
    caption = "Top 5 UP regulated Genes in Discrepancy Cells. filter(p_val_adj < 0.05)",
    digits = c(0, 3, 3, 3, 2, 2), 
    col.names = c("Gene", "Avg Log2FC", "pct.1 in Disagree", "pct.2 in Agree", "p-value", "Adj p-value")
  ) %>%
  add_header_above(c(" " = 1, "Comparison: ident.1 = Disagree vs ident.2 = Agree" = 5)) %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"), full_width = FALSE)
```


```{r}
discrepancy_up_DEGs <- FindMarkers(
  query,
  ident.1 = "Disagree",
  ident.2 = "Agree",
  only.pos = TRUE,         
  min.pct = 0.10,          
  logfc.threshold = 0.25)

top5_UP_no_filter <- discrepancy_up_DEGs %>%
  tibble::rownames_to_column(var = "Gene") %>%
  arrange(desc(avg_log2FC)) %>%
  head(5) %>%
  select(Gene, avg_log2FC, pct.1, pct.2, p_val, p_val_adj)

top5_UP_no_filter %>%
  kable(
    caption = "Top 5 UP regulated Genes in Discrepancy Cells. NO filter(p_val_adj < 0.05)",
    digits = c(0, 3, 3, 3, 2, 2), 
    col.names = c("Gene", "Avg Log2FC", "pct.1 in Disagree", "pct.2 in Agree", "p-value", "Adj p-value")
  ) %>%
  add_header_above(c(" " = 1, "Comparison: ident.1 = Disagree vs ident.2 = Agree" = 5)) %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"), full_width = FALSE)
```


```{r}
discrepancy_down_DEGs <- FindMarkers(
  query,
  ident.1 = "Disagree",
  ident.2 = "Agree",
  only.pos = FALSE,         
  min.pct = 0.10,          
  logfc.threshold = 0.25)

top5_DOWN_filter <- discrepancy_down_DEGs %>%
  tibble::rownames_to_column(var = "Gene") %>%
  filter(avg_log2FC < 0 & p_val_adj < 0.05) %>%  
  arrange(avg_log2FC) %>%                          
  head(5) %>%
  select(Gene, avg_log2FC, pct.1, pct.2, p_val, p_val_adj)


top5_DOWN_filter %>%
  kable(
    caption = "Top 5 DOWN regulated Genes in Discrepancy Cells with filter(p_val_adj < 0.05)",
    digits = c(0, 3, 3, 3, 2, 2),
    col.names = c("Gene", "Avg Log2FC", "pct.1 in Disagree", "pct.2 in Agree", "p-value", "Adj p-value")
  ) %>%
  add_header_above(c(" " = 1, "Comparison: ident.1 = Disagree vs ident.2 = Agree" = 5)) %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"), full_width = FALSE)

```


```{r}

discrepancy_down_DEGs <- FindMarkers(
  query,
  ident.1 = "Disagree",
  ident.2 = "Agree",
  only.pos = FALSE,         
  min.pct = 0.10,          
  logfc.threshold = 0.25)

top5_DOWN_no_filter <- discrepancy_down_DEGs %>%
  tibble::rownames_to_column(var = "Gene") %>%        
  arrange(avg_log2FC) %>%                              
  head(5) %>%
  select(Gene, avg_log2FC, pct.1, pct.2, p_val, p_val_adj)


top5_DOWN_no_filter %>%
  kable(
    caption = "Top 5 DOWN regulated Genes in Discrepancy Cells with NO filter(p_val_adj < 0.05)",
    digits = c(0, 3, 3, 3, 2, 2),
    col.names = c("Gene", "Avg Log2FC", "pct.1 in Disagree", "pct.2 in Agree", "p-value", "Adj p-value")
  ) %>%
  add_header_above(c(" " = 1, "Comparison: ident.1 = Disagree vs ident.2 = Agree" = 5)) %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"), full_width = FALSE)
```


```{r}

top5_filtered <- discrepancy_down_DEGs %>%
  tibble::rownames_to_column(var = "Gene") %>%
  filter(p_val_adj < 0.05 & avg_log2FC < 0) %>%
  arrange(avg_log2FC) %>%
  head(5) %>%
  select(Gene, avg_log2FC, pct.1, pct.2, p_val, p_val_adj)

top5_unfiltered <- discrepancy_down_DEGs %>%
  tibble::rownames_to_column(var = "Gene") %>%
  filter(avg_log2FC < 0) %>%
  arrange(avg_log2FC) %>%
  head(5) %>%
  select(Gene, avg_log2FC, pct.1, pct.2, p_val, p_val_adj)

t1 <- top5_filtered %>%
  kable(
    caption = "<b>Top Fold-Change DOWN regulated Genes (filter: p_val_adj < 0.05)</b>",
    digits = c(0, 3, 3, 3, 2, 2),
    col.names = c("Gene", "Avg Log2FC", "pct.1 (Disagree)", "pct.2 (Agree)", "p-value", "Adj p-value")
  ) %>%
  add_header_above(c(" " = 1, "Comparison: ident.1 = Disagree vs ident.2 = Agree" = 5)) %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"), full_width = FALSE)

t2 <- top5_unfiltered %>%
  kable(
    caption = "<b>Top Fold-Change DOWN regulated Genes (NO filter)</b>",
    digits = c(0, 3, 3, 3, 2, 2),
    col.names = c("Gene", "Avg Log2FC", "pct.1 (Disagree)", "pct.2 (Agree)", "p-value", "Adj p-value")
  ) %>%
  add_header_above(c(" " = 1, "Comparison: ident.1 = Disagree vs ident.2 = Agree" = 5)) %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"), full_width = FALSE)

kables(list(t1, t2))
```


```{r}

top_genes <- rownames(top10_discrepancy_genes)

top_4_discrepancy_genes <- FeaturePlot(
  query,
  features = top_genes[1:min(4, length(top_genes))],
  reduction = "umap",
  cols = c("grey", "blue"),
  ncol = 2
) & theme_minimal()

top_4_discrepancy_genes

ggsave( "/Users/kiyounghan/Desktop/top_4_discrepancy_genes.png", plot = top_4_discrepancy_genes,width  = 20 ,height = 15,dpi = 400)
```


```{r}

vlnplot_agree_disagree <- VlnPlot(
  query,
  features = top_genes[1:min(4, length(top_genes))],
  group.by = "model_agreement",
  cols = c("Agree" = "grey", "Disagree" = "red"),
  pt.size = 0.0,
  ncol = 2
) & 
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) &
  theme(legend.position = "none")

vlnplot_agree_disagree 

ggsave( "/Users/kiyounghan/Desktop/vlnplot_agree_disagree.png", plot = vlnplot_agree_disagree, width  = 20 ,height = 15,dpi = 400)

```

print("========================================================================================================")
print(" 11. Cross Evaluation: Random Forest vs Seurat...COMPLETED                                              ")
print("========================================================================================================")



print("===================================================")
print(" 12. EXPRESSION HEATMAP ... STARTING               ")
print("===================================================")

```{r}

Idents(query) <- "rf_model_Labels_clean"
cat("Finding markers across Random Forest clusters...\n")
rf_all_markers <- FindAllMarkers(
  query,
  only.pos        = TRUE,
  min.pct         = 0.25,
  logfc.threshold = 0.25 )

top10_rf_markers <- rf_all_markers %>%
  group_by(cluster) %>%
  filter(p_val_adj < 0.05) %>%
  slice_max(n = 10, order_by = avg_log2FC)

print(top10_rf_markers %>% select(cluster, gene, avg_log2FC, p_val_adj))
write.csv(top10_rf_markers, file = "top10_rf_markers.csv", row.names = FALSE)
top_gene <- unique(top10_rf_markers$gene)
query <- ScaleData(query, features = top_gene)

```


```{r}
DefaultAssay(query) <- "RNA"

avg_exp <- AverageExpression(
  query,
  features = top_gene,
  group.by = "rf_model_Labels_clean",
  slot     = "data"
)$RNA
```


```{r}
scaled_avg_exp <- t(scale(t(avg_exp)))

temp_hp <- pheatmap(scaled_avg_exp, cluster_cols = TRUE, silent = TRUE)
clustered_col_order <- colnames(scaled_avg_exp)[temp_hp$tree_col$order]

top10_rf_markers$cluster <- factor(top10_rf_markers$cluster, levels = clustered_col_order)

genes_ordered <- top10_rf_markers %>%
  filter(!is.na(cluster)) %>%
  arrange(cluster) %>%
  pull(gene) %>%
  unique()
```


```{r}
scaled_avg_sub <- scaled_avg_exp[genes_ordered, clustered_col_order]

Expression_Heatmap <- pheatmap(
  scaled_avg_sub,
  cluster_rows = FALSE,
  cluster_cols = TRUE,
  show_colnames = TRUE,
  show_rownames = TRUE,
  border_color = NA,
  fontsize_row = 4.5,
  fontsize_col = 11,
  angle_col    = 45,
  width        = 9,
  height       = 11,
  color = colorRampPalette(c("blue", "white", "red"))(100)
  )

Expression_Heatmap

ggsave( "/Users/kiyounghan/Desktop/Expression_Heatmap.png", plot = Expression_Heatmap,width  = 20 ,height = 15,dpi = 400)

```

print("===================================================")
print(" 12. EXPRESSION HEATMAP ... COMPLETED              ")
print("===================================================")







print("=============================================================================================================")
print(" I I I. SUPERVISED (Reference: ref ("pbmcsca") & Query: query ("pbmc3k")) by naive Random Forest...COMPLETED ")
print("=============================================================================================================")









print("============================================================================================")
print(" I V. Pseudo-bulk analysis for single-cell RNA-Seq data (Data = "pbmcsca" as pbmcsca)       ")
print("============================================================================================")


```{r}
rm(list = ls())
```


```{r setup, include=FALSE}
pacman::p_load(apeglm, BiocManager, caret, clue, DESeq2, doParallel, dplyr, e1071, EnhancedVolcano, fgsea, foreach, ggplot2, gridExtra, gt, irr, kableExtra, knitr, Matrix, mclust, msigdbr, parallel, patchwork, pheatmap, randomForest, ranger, Seurat, SeuratData, tidyverse )
```


```{r}
pbmcsca <- LoadData("pbmcsca")
pbmcsca <- UpdateSeuratObject(pbmcsca)
```


### single Cell


```{r}
head(pbmcsca@meta.data,3)
```


```{r}
possible_names <- c("cell_type", "celltype", "CellType", "cell.type", "labels")
found_name <- intersect(possible_names, colnames(pbmcsca@meta.data))

if (length(found_name) > 0) {
  pbmcsca$CellType <- pbmcsca@meta.data[[found_name[1]]]
  message("Success!!! Found and copied cell labels from metadata column: ==>>'", found_name[1], "'")
} else {
  print(colnames(pbmcsca@meta.data))
  stop("Could not find cell type metadata. Look at the printed column names above and choose the right one!")
}
```


```{r}
factor(sort(as.numeric(unique(pbmcsca@meta.data$Cluster))))
```


```{r}
unique(pbmcsca@meta.data$orig.ident)
```


```{r}
unique(pbmcsca@meta.data$Experiment)
```


```{r}
unique(pbmcsca@meta.data$Method)
```


```{r}
unique(pbmcsca@meta.data$CellType)
```


```{r}
singlecell_counts <- GetAssayData(pbmcsca, assay = "RNA",  layer = "counts")[1:3,]
singlecell_counts
```



```{r}
singlecell_DEGs <- FindMarkers(pbmcsca , 
                            group.by = "CellType",   
                       ident.1 = "Cytotoxic T cell" ,
                    ident.2 = "CD14+ monocyte",
                       test.use = "wilcox")

print(
  nrow(
    singlecell_DEGs[singlecell_DEGs$p_val_adj < 0.05,]
      )  )

```


```{r}
head(singlecell_DEGs)
```








############# Pseudobulk #############################


```{r}
pseudobulk_matrix <- AggregateExpression(
           pbmcsca, 
         assays = "RNA", 
        return.seurat = FALSE, 
  group.by = c("CellType", "orig.ident"), 
  slot = "counts"
)$RNA


dim(pseudobulk_matrix)
head(pseudobulk_matrix[1:3, 1:7])
```


```{r}
class(pseudobulk_matrix)
```


```{r}
selected_cells <- colnames(pbmcsca)[
  pbmcsca$CellType == "B cell" & pbmcsca$orig.ident == "pbmc1"
]

sum(singlecell_counts["TSPAN6", selected_cells])

pseudobulk_matrix["TSPAN6","B cell_pbmc1"]
```



```{r}
df_info <- data.frame( 
  df_id = colnames(pseudobulk_matrix) ,
   row.names = colnames(pseudobulk_matrix) )

#  extract Cell Type
df_info$CellType <- sub("_[^_]+$", "", df_info$df_id)

# extract orig.ident
df_info$orig.ident <- sub(".*_", "", df_info$df_id)  
rownames(df_info) <- df_info$df_id

# check for alignment
identical(colnames(pseudobulk_matrix), rownames(df_info))

```


```{r}
head(df_info)
```


```{r}

target_celltypes <- c("Cytotoxic T cell", "CD4+ T cell")
meta_sub <- df_info %>% filter(CellType %in% target_celltypes)
meta_sub
```



```{r}
common_samples <- intersect(colnames(pseudobulk_matrix), rownames(meta_sub))
common_samples
class(common_samples)

```


```{r}
meta_sub <- meta_sub[common_samples, ]
meta_sub
```


```{r}
counts_sub <- pseudobulk_matrix[, rownames(meta_sub)]
head(counts_sub) 
```


```{r}
meta_sub$CellType   <- factor(meta_sub$CellType, levels = target_celltypes)
meta_sub$orig.ident <- factor(meta_sub$orig.ident)
```


```{r}
meta_sub$CellType
meta_sub$orig.ident
class(c(meta_sub$CellType,meta_sub$orig.ident))
```


```{r}
DE_seq_object <- DESeqDataSetFromMatrix(
  countData = as.matrix(counts_sub),
  colData   = meta_sub,
  design    =  ~ orig.ident + CellType
)
```


```{r}
var_stab_transf <- vst(DE_seq_object, blind = FALSE)
head(assay(var_stab_transf))
```


```{r}
2^(13.066656 - 10.205740)
2^(13.587725 - 8.329803) 
```


```{r}
PCA_var_stab_transf <- plotPCA(var_stab_transf, intgroup = c("orig.ident", "CellType"))
PCA_var_stab_transf
```


```{r}
DE_seq_object <- DE_seq_object[rowSums(counts(DE_seq_object)) >= 10, ]

DE_seq_object <- DESeq(DE_seq_object)
```


```{r}
DE_seq_object$sizeFactor
```


```{r}
head(results(DE_seq_object))
```


```{r}
DE_seq_object_result <- as.data.frame(results(DE_seq_object, contrast = c("CellType", "Cytotoxic T cell", "CD4+ T cell"))) %>%
  rownames_to_column("gene") %>% 
  arrange(padj)

head(DE_seq_object_result)
```



# lfcShrink()
```{r}
resultsNames(DE_seq_object)

apeglm_result <- lfcShrink(DE_seq_object, coef = resultsNames(DE_seq_object)[3], type = "apeglm")
```


```{r}
apeglm_df <- as.data.frame(apeglm_result) %>%
  select(-any_of("gene")) %>%  
  rownames_to_column("gene") %>%
  arrange(padj)

head(apeglm_df)
```

```{r}
DE_seq_object
```


```{r}
head(counts(DE_seq_object))
```


```{r}
colData(DE_seq_object)
```


```{r}
apeglm_result <- as.data.frame(apeglm_result)
apeglm_result$gene <- rownames(apeglm_result)

significant_diff_exp_genes <- subset(apeglm_result, padj < 0.05 & abs(log2FoldChange) >= 1.0)

significant_diff_exp_genes <- significant_diff_exp_genes[order(significant_diff_exp_genes$padj), ]

top20_significant_diff_exp_genes <- head(significant_diff_exp_genes,20)
top20_significant_diff_exp_genes
```


```{r}
gene_names <- top20_significant_diff_exp_genes$gene 
top20_genes <- assay(var_stab_transf)[gene_names, ]

top20_genes <- as.matrix(top20_genes)

all(apply(top20_genes, 2, is.numeric))
```


```{r}
annotation_color <- data.frame(
  CellType = vsd$CellType,
  Donor = vsd$orig.ident,
  row.names = colnames(vsd)
)

pheatmap(top20_genes, 
         scale = "row", 
         annotation_col = annotation_color,
         show_colnames = TRUE,
         treeheight_row = 30,
         treeheight_col = 30,
         main = "Top 20 DEGs (Cytotoxic T Cells vs CD4+)"
         )
```
