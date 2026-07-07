print("=======================================================================================") 
print("  1 a. Libraries... STARTING                                                           ")  
print("=======================================================================================")

library(caret)
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
if(!require(e1071)) install.packages("e1071", repos="http://cran.us.r-project.org")
library(e1071)

print("=======================================================================================")
print("  1 b. Libraries... COMPLETED                                                          ")
print("=======================================================================================")



print("========================================================================================")
print(" 1 b. Upload spare data.... STARTING                                                    ") 
print("========================================================================================")

data_path <- "/Users/kiyounghan/Desktop/pbmc2700/filtered_gene_bc_matrices/hg19"
pbmc.data <- Read10X(data.dir = data_path)
pbmc_2700 <- CreateSeuratObject(counts = pbmc.data, project = "pbmc2700", min.cells = 3, min.features = 200)

print("=======================================================================================")
print(" 1 b. Upload spare data... COMPLETED                                                   ")
print("=======================================================================================")




print("=======================================================================================")
print("  2. Data cleaning... STARTING                                                         ")
print("=======================================================================================")

pbmc_2700[["percent.mt"]] <- PercentageFeatureSet(pbmc_2700, pattern = "^MT-")
pbmc_2700 <- subset(pbmc_2700, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc_2700 <- NormalizeData(pbmc_2700, verbose = FALSE)
pbmc_2700 <- FindVariableFeatures(pbmc_2700, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
pbmc_2700 <- ScaleData(pbmc_2700, verbose = FALSE)
pbmc_2700 <- RunPCA(pbmc_2700, features = VariableFeatures(object = pbmc_2700), verbose = FALSE)

print("=======================================================================================")
print(" 2. Data cleaning... COMPLETED                                                         ")
print("=======================================================================================")




print("=======================================================================================")
print(" 3. Clustering & Dimension Reduction... Starting                                       ")
print("=======================================================================================")

# --- METHOD A. K-Means --- 

query_pca <- Embeddings(pbmc_2700, reduction = "pca")[, 1:10]
set.seed(123)
message("Evaluating optimal cluster count (k) using the Elbow Method...")
wcss <- vector()
for (k in 1:10) {
  km <- kmeans(query_pca, centers = k, nstart = 10)
  wcss[k] <- km$tot.withinss
}

plot(
  1:10, wcss,
  type = "b",
  pch  = 19,
  col  = "darkblue",
  xlab = "Number of Clusters (k)",
  ylab = "Within-Cluster Sum of Squares (WCSS)",
  main = "Elbow Method for Optimal K-means"
)

# 2nd Derivative - concavity
wcss_derivatives <- diff(diff(wcss))
optimal_k <- which.max(wcss_derivatives) + 1

message(paste("Mathematical elbow detected at k =", optimal_k))
kmeans_output <- kmeans(query_pca, centers = optimal_k, nstart = 10)  
pbmc_2700$KMeans_Clusters <- as.factor(kmeans_output$cluster)
  
print("=========================================================================================================")
print("Math (Elbow method) shows 2 clusters: Lymphoid(T,B,NK Cells) vs Myeloid(Monocytes,Dendritic Cells)       ")
print("Biology dictates 4 clusters: Lymphoid should be split into B-Cells MS4A1, T-Cells CD3D, NK-Cells GNL     ")
print("=========================================================================================================")


# --- METHOD B: Hierarchical Clustering  ---

distance_matrix <- dist(query_pca)
hierarchical_tree <- hclust(distance_matrix, method = "ward.D2")
pbmc_2700$Hierarchical_Clusters <- as.factor(cutree(hierarchical_tree,optimal_k ))     
message(paste("Hierarchical tree cut into", optimal_k, "clusters."))   

message("Evaluating k-means vs hierarchical cluster alignment:")
print(table(KMeans = pbmc_2700$KMeans_Clusters, Hierarchical = pbmc_2700$Hierarchical_Clusters))

aa <- sum( 
  ( table(KMeans = pbmc_2700$KMeans_Clusters, Hierarchical = pbmc_2700$Hierarchical_Clusters)[1,1] + 
    table(KMeans = pbmc_2700$KMeans_Clusters, Hierarchical = pbmc_2700$Hierarchical_Clusters)[2,2] ) / 
    sum(table(KMeans = pbmc_2700$KMeans_Clusters, Hierarchical = pbmc_2700$Hierarchical_Clusters))
) * 100

message(paste("K-Means and Hierarchical clustering share" ,aa,"% alignment on this dataset."))
library(raster)
kappacoefficient <- kappa(table(KMeans = pbmc_2700$KMeans_Clusters, Hierarchical = pbmc_2700$Hierarchical_Clusters))
kappacoefficient





print("===========================================================================================")
print("     Unconstrained Louvain... Starting                                                     ")
print("===========================================================================================")

# --- METHOD C: Louvain ---

pbmc_2700 <- FindNeighbors(pbmc_2700, dims = 1:10)
pbmc_2700 <- FindClusters(pbmc_2700, resolution = 0.5, verbose = FALSE)

# Extract
actual_clusters    <- sort(as.numeric(as.character(unique(Idents(pbmc_2700)))))
num_clusters_found <- length(actual_clusters)

message(paste("Louvain converged on", num_clusters_found, "distinct clusters."))
message(paste("Cluster IDs discovered:", paste(actual_clusters, collapse = ", ")))

# Dictionary
pbmc_dictionary <- c(
  "0" = "Naive CD4+ T Cells",
  "1" = "CD14+ Monocytes",
  "2" = "Memory CD4+ T Cells",
  "3" = "B Cells",
  "4" = "CD8+ T Cells",
  "5" = "FCGR3A+ Monocytes",
  "6" = "NK Cells",
  "7" = "Dendritic Cells",
  "8" = "Platelets"
)

# Subset dictionary to equal what I found
cell_identities <- pbmc_dictionary[as.character(actual_clusters)]

if (any
        (is.na
                (cell_identities)
)
)
 {
  missing_indices <- which(is.na(cell_identities))
  for (i in missing_indices) {
    cluster_num <- actual_clusters[i]
    cell_identities[i] <- paste("Unclassified_Cluster_", cluster_num, sep="")
  }
}
  
# Rename INSIDE the Seurat
pbmc_2700 <- RenameIdents(pbmc_2700, cell_identities)
  
# Save names to a metadata column for later use
pbmc_2700$Louvain_Biological_Names <- Idents(pbmc_2700)

#message("Biological identities mapped.")
#print(table(pbmc_2700$Louvain_Biological_Names))

print("========================================================================================")
print(" 3. Clustering & Dimension Reduction.... COMPLETED                                      ")
print("========================================================================================")





print("========================================================================================")
print("  4. Visualization.... STARTING                                                         ")
print("========================================================================================")

pbmc_2700 <- RunUMAP(pbmc_2700, dims = 1:10, verbose = FALSE)
pbmc_2700 <- RunTSNE(pbmc_2700, dims = 1:10, verbose = FALSE)

p1 <- DimPlot(pbmc_2700, reduction = "pca", group.by = "KMeans_Clusters", pt.size = 0.3) +  ggtitle("1: K-Means")
p2 <- DimPlot(pbmc_2700, reduction = "pca", group.by = "Hierarchical_Clusters", pt.size = 0.3) + ggtitle("2: Hierarchical")
p3 <- DimPlot(pbmc_2700, reduction = "tsne", pt.size = 0.3) + ggtitle("3: t-SNE")
p4 <- DimPlot(pbmc_2700, reduction = "umap", pt.size = 0.3) + ggtitle("4: UMAP")
        
message("4 figures panel")
library(patchwork)
final_plot <- (p1 + p2) / (p3 + p4)

ggsave( "/Users/kiyounghan/Desktop/clustering_evolution_progression.png",plot   = final_plot,width  = 12,height = 10,dpi = 300)

print("================================================================================================================")
print(" 4. Visualization COMPLETED. Check DESKTOP for /Users/kiyounghan/Desktop/clustering_evolution_progression.png   ")
print("================================================================================================================")





print("====================================================================================")
print("  5.  Lineage Verification Plots (Canonical Marker Expression) STARTING             ")   
print("====================================================================================") 

marker_plot <- FeaturePlot(
        pbmc_2700                                     ,
        features = c("MS4A1", "CD14", "CD8A", "GNLY") ,
        ncol = 2                                      ,
        reduction = "umap"
)
 
marker_plot

ggsave(
        "/Users/kiyounghan/Desktop/canonical_markers_verification.png" ,
         plot = marker_plot                                  ,
         width = 10                                          ,
         height = 8                                          ,
         dpi = 300
)

print("====================================================================================")
print(" 5.Lineage Verification Plots(Canonical Marker Expression) COMPLETED                ")   
print("====================================================================================")




print("============================================================================================")
print(" I I.    SUPERVISED (Reference: pbmc16 & Query: pbmc_2700) SVM (KERNEL = RBF) Starting      ")
print("============================================================================================")


print("========================================================================================")
print(" 1. Upload spare data & prep.... STARTING                                               ") 
print("========================================================================================")

# --- QUERY DATA ---

query_path <- "/Users/kiyounghan/Desktop/pbmc2700/filtered_gene_bc_matrices/hg19"
query_data <- Read10X(data.dir = query_path)
query <- CreateSeuratObject( counts = query_data, min.cells = 3, min.features = 200 )


# --- REFERENCE DATA ---

library(devtools)
#devtools::install_github("satijalab/seurat-data")
library(Seurat)
library(SeuratData)

#readRDS("/Users/kiyounghan/Desktop/manifest.rds")
#AvailableData()

pbmcsca_data <- SeuratData::LoadData("pbmcsca")
pbmcsca_data <- subset(pbmcsca_data, subset = Method == "10x Chromium (v2)")

# (1 )SEURAT v5 SAFE: Pull the count layer directly using v5 bracket syntax
# works on 'Assay5' objects and avoids all slot errors!
# ref_counts is just a raw numeric matrix
ref_counts <- pbmcsca_data[["RNA"]]$counts


# 2. CANNOT!!!!! do metadata assignment before CreateSeuratObject()
#    Seurat Object is designed explicitly to keep data synchronized. 
#    @meta.data slot (dedicated data frame) where every row maps to a column(cell) in expression matrix.
ref <- CreateSeuratObject(counts = ref_counts)



# 3. METADATA ASSIGNMENT: Auto-detect the cell type column name
possible_names <- c("cell_type", "celltype", "CellType", "cell.type", "labels")
found_name <- intersect(possible_names, colnames(pbmcsca_data@meta.data))

if (length(found_name) > 0) {
  # Grab the first matching column name found in the metadata matrix
  ref$CellType <- pbmcsca_data@meta.data[[found_name[1]]]
  message("Success! Found and copied cell labels from metadata column: '", found_name[1], "'")
} else {
  # Fallback: list all columns so can see what it's called
  print(colnames(pbmcsca_data@meta.data))
  stop("Could not find cell type metadata. Look at the printed column names above and choose the right one!")
}

print("=======================================================================================")
print(" 1. Upload spare data & prep... COMPLETED                                              ")
print("=======================================================================================")






print("====================================================================================")
print(" 2. VERIFICATION : Intersection genes...starting                                    ")
print("====================================================================================")

shared_genes <- intersect(rownames(ref), rownames(query))
print(paste("Number of shared genes found:", length(shared_genes)))
print( head(shared_genes , 7 )  )

print("====================================================================================")
print(" 2. VERIFICATION : Intersection genes...COMPLETED                                   ")
print("====================================================================================")





print("====================================================================================")
print(" 3. NOTE for future use...starting                                                  ")
print("====================================================================================")

print(colnames(pbmcsca_data@meta.data))
print(sort(as.numeric(unique(pbmcsca_data@meta.data$Cluster))))

print("====================================================================================")
print(" 3. NOTE for future use...COMPLETED                                                 ")
print("====================================================================================")






print("====================================================================================")
print(" 4. Seurat Objects, Normalize, FindVariableFeatures ...starting                     ")
print("====================================================================================")

ref <- NormalizeData(ref, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
ref <- FindVariableFeatures(ref, selection.method = "vst", nfeatures = 2000, verbose = FALSE)

# Dont do CellType bc we want to train SVM on reference and assign to query so CANT have name assigned.
query <- CreateSeuratObject(counts = query_data)
query <- NormalizeData(query, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
query <- FindVariableFeatures(query, selection.method = "vst", nfeatures = 2000, verbose = FALSE)

print("====================================================================================")
print(" 4. Seurat Objects, Normalize, FindVariableFeatures...COMPLETED                     ")
print("====================================================================================")





print("====================================================================================")
print(" 5. Intersection genes  ...starting                                                 ")
print("====================================================================================")

ref_genes <- VariableFeatures(ref)
query_genes <- VariableFeatures(query)

# hvgs is commonly used for "Highly Variable Genes". 
# genes that are statistically filtered b/c their variance across cells is significantly higher than background noise.
# high-signal, high-variance biological marker genes driving PCA and SVM

shared_hvgs <- intersect(ref_genes, query_genes)
final_features <- head(shared_hvgs, 2000)

common_features <- intersect(final_features, rownames(query[["RNA"]]))

print("====================================================================================")
print(" 5. Intersection genes...COMPLETED                                                  ")
print("====================================================================================")





print("====================================================================================")
print(" 6. Data management for Seurat, caret, & RF compatibility...starting                ")
print("====================================================================================")

# Seurat(rows=genes & columns=cells) RF(rows=cells & columns=genes) 

# select common_features= hgvs & drop uninformative genes
# sparse matrix into R matrix. MAYBE INEFFICIENT METHOD but the only way I know.
ref_counts <- t(as.matrix(GetAssayData(ref, layer = "data")[common_features, ]))
ref_df     <- as.data.frame(ref_counts)
ref_df$Label <- as.factor(ref$CellType)

# Query
query_counts <- t(as.matrix(GetAssayData(query, layer = "data")[common_features, ]))
query_features <- as.data.frame(query_counts)

# Fix column names so caret doesn't crash on gene symbols (like '-' or '.')
# MAJOR pain !!! BE CAREFUL IN THE FUTURE
colnames(ref_df)[1:length(common_features)] <- make.names(colnames(ref_df)[1:length(common_features)])
colnames(query_features)[1:length(common_features)] <- make.names(colnames(query_features)[1:length(common_features)])

print("====================================================================================")
print(" 6. Data management for Seurat, caret, & RF compatibility...COMPLETED               ")
print("====================================================================================")





print("====================================================================================")
print(" 7.  80 % Train / 20 % Test STRATIFIED!! Sampling...starting                        ")
print("====================================================================================")

library(dplyr)
set.seed(123)

train_idx       <- caret::createDataPartition(ref_df$Label, p = 0.80, list = FALSE)
train_pool_full <- ref_df[train_idx, ]
test_pool       <- ref_df[-train_idx, ]

# max 300 cells per type
set.seed(123)
train_pool <- train_pool_full %>%
  dplyr::group_by(Label) %>%
  dplyr::slice_sample(n = 300, replace = FALSE) %>%
  dplyr::ungroup()

train_pool <- as.data.frame(train_pool)

# Ensure clean syntax rules for the class labels
train_pool$Label         <- as.factor(train_pool$Label)
test_pool$Label          <- as.factor(test_pool$Label)
levels(train_pool$Label) <- make.names(levels(train_pool$Label))
levels(test_pool$Label)  <- make.names(levels(test_pool$Label))

print("====================================================================================")
print(" 7.  80 % Train / 20 % Test STRATIFIED!! Sampling...COMPLETED                       ")
print("====================================================================================")









print("====================================================================================")
print(" 8. Parallel Processors & 10-fold CV...starting                                     ")            
print("====================================================================================")

library(parallel)
library(doParallel)
library(kernlab)

cores_use <- 10
cl <- makeCluster(cores_use)
registerDoParallel(cl)

rf_grid <- expand.grid(
  mtry = c(2, 5, 15, 30, 50, 100)
)

# "allowParallel = TRUE" => caret to do 10-fold using parallel computing

cv_control <- caret::trainControl(
  method          = "cv"     ,
  number          = 5       ,
  savePredictions = "final"  ,
  allowParallel   = TRUE     ,
  classProbs      = TRUE
 )

# metric = "Kappa" optimizes it and Accuracy
# importance = TRUE calculates which genes matter most
set.seed(123)
rf_model <- caret::train(
  Label ~ .,
  data      = train_pool,
  method    = "rf",          # Switched to Random Forest!
  trControl = cv_control,
  tuneGrid  = rf_grid,
  ntree     = 250,
  metric    = "Kappa",
  importance = TRUE
)

#  IMPORTANT!!! : Shut down the cluster immediately after training ends
stopCluster(cl)
registerDoSEQ()

message("Training complete. CPU back to single core.")
print(rf_model$bestTune)

print("====================================================================================")
print(" 8. Parallel Processsors & 10-fold CV...COMPLETED                                   ")
print("====================================================================================")







print("====================================================================================")
print(" 9. RF Model Accuracy on Test Pool...starting                                       ")
print("====================================================================================")

# Isolate the feature columns from the reference's test pool and leave out Label column.
test_features    <- test_pool[, 1:length(common_features)]
test_predictions <- predict(rf_model, newdata = test_features)
conf_matrix      <- confusionMatrix(test_predictions, test_pool$Label)

# Kappa Statistics is between -1 & 1 with 0 meaning noise but the statistic is not accepted all the time. 
print(paste("True Validation Accuracy:", round(conf_matrix$overall["Accuracy"], 4)))
print(paste("True Kappa Statistic:", round(conf_matrix$overall["Kappa"], 4)))

print("====================================================================================")
print(" 9. RF Model Accuracy on Test Pool...COMPLETED                                      ")
print("====================================================================================")




print("====================================================================================")
print(" 10. Query Projection & Prediction...starting                                       ")
print("====================================================================================")

predicted_classes <- predict(rf_model, newdata = query_features)

# IMPORTANT SECTION !!! PREVENTS BAD predictions
final_predictions <- as.character(predicted_classes)
final_predictions <- gsub("\\.\\.", "+ ", final_predictions)
final_predictions <- gsub("\\.", " ", final_predictions)

query$rf_model_Labels <- final_predictions

#print(" *** FINAL RANDOM FOREST CLASSIFICATIONS *** ")
print(table(query$rf_model_Labels))

print("====================================================================================")
print(" 10. Query Projection & Prediction...COMPLETED                                      ")
print("====================================================================================")









print("===============================================================================================")
print(" 11. Cross Evaluation: Random Forest vs Louvain    ...starting                                 ")
print("===============================================================================================")

q_cells   <- unname(as.character(rownames(query@meta.data)))
r_cells   <- unname(as.character(rownames(pbmc_2700@meta.data)))
r_labels  <- unname(as.character(pbmc_2700@meta.data$Louvain_Biological_Names))

# 2. Match positions
indices <- match(q_cells, r_cells)
final_labels <- r_labels[indices]

# 3. CRITICAL: Force fill ALL missing values right now
final_labels[is.na(final_labels)] <- "Other"
final_labels[final_labels == ""]  <- "Other" # Handle empty strings just in case

#print(table(final_labels, useNA = "always"))

# 4. Lock it into the Seurat metadata frame
query@meta.data[["Louvain_Biological_Names"]] <- as.factor(final_labels)
query@meta.data[["rf_model_Labels"]]          <- as.factor(unname(as.character(query@meta.data[["rf_model_Labels"]])))


contingency_table <- table(
  Louvain = query@meta.data[["Louvain_Biological_Names"]], 
  RF      = query@meta.data[["rf_model_Labels"]]
)

print(contingency_table)

# Row sum to 100 %
print(round(prop.table(contingency_table, margin = 1) * 100, 1))

# Save to CSV file
write.csv(as.data.frame.matrix(contingency_table), file = "RF_vs_Louvain_contingency_table.csv")

print("==================================================================================================================")
print(" 11. Cross Evaluation: Random Forest vs Louvain...COMPLETED                                                      ")
print("==================================================================================================================")







print("=========================================================================================")
print(" 12. Visualization...starting                                                            ")
print("=========================================================================================")

# ARGHHHHHHH cut off any spaces in front and back of barcodes. WHY  ARE THEY HERE 
ref_barcodes   <- trimws(rownames(pbmc_2700@meta.data))
query_barcodes <- trimws(rownames(query@meta.data))

# 2. Extract reference labels as characters
louvain_labels <- as.character(pbmc_2700@meta.data$Louvain_Biological_Names)

# 3. Match positions
matching_indices <- match(query_barcodes, ref_barcodes)

# 4. Extract matched values into a temporary character vector (no factor rules yet)
temp_labels <- louvain_labels[matching_indices]

# 5. FIXED: The assignment arrow actually overwrites the NAs with "Other"
temp_labels[is.na(temp_labels)] <- "Other"

# 6. Strip vector names to prevent structural metadata conflicts
temp_labels <- unname(temp_labels)

# 7. Convert to factors and force directly into the metadata frame
query@meta.data[["Louvain_Biological_Names"]] <- as.factor(temp_labels)
query@meta.data[["rf_model_Labels"]]          <- as.factor(unname(as.character(query@meta.data[["rf_model_Labels"]])))


# Normalize, FindVariableFeatures, ScaleData, PCA, UMAP on query  
query <- NormalizeData(query, normalization.method = "LogNormalize", scale.factor = 10000)
query <- FindVariableFeatures(query, selection.method = "vst", nfeatures = 2000)
query <- ScaleData(query, features = VariableFeatures(query))
query <- RunPCA(query, features = VariableFeatures(query), verbose = FALSE)
query <- RunUMAP(query, dims = 1:30, reduction = "pca", verbose = FALSE)


Idents(query) <- "Louvain_Biological_Names"
p1 <- DimPlot(query, reduction = "umap", label = TRUE, repel = TRUE) + 
  ggtitle("Unsupervised Louvain Clusters") + 
  theme(legend.position = "none")


Idents(query) <- "rf_model_Labels"
p2 <- DimPlot(query, reduction = "umap", label = TRUE, repel = TRUE) + 
  ggtitle("Supervised Random Forest Labels") + 
  theme(legend.position = "none")


# side-by-side in color
options(repr.plot.width = 12, repr.plot.height = 6)
final_plot = p1 + p2
ggsave( "/Users/kiyounghan/Desktop/Louvain_Random_Forest.png", plot = final_plot,width  = 12,height = 10,dpi = 300) 

print("=========================================================================================")
print(" 12. Visualization...COMPLETED                                                           ")
print("=========================================================================================")





print("==============================================================================")
print(" 13. Downstream Diagnostics...starting                                        ")
print("==============================================================================")

# --- Top misclassifications - Discrepancy Analysis ---

louvain_text <- as.character(query$Louvain_Biological_Names)
rf_text      <- as.character(query$rf_model_Labels)

# 2. Optional: Standardize common text naming mismatches if needed
# (e.g., if one says "B Cells" and the other says "B cell")
louvain_text <- tolower(gsub("s$", "", louvain_text)) # lowercase and remove trailing 's'
rf_text      <- tolower(gsub("s$", "", rf_text))      # lowercase and remove trailing 's'

# 3. Create the logical mismatch column (TRUE if they disagree, FALSE if they agree)
query$label_mismatch <- louvain_text != rf_text

# See on map. 
message("Plotting discrepancies on UMAP...")

# Highlight the mismatched cells (Mismatches in dark red, agreements in light gray)
p1 <- DimPlot(query, group.by = "label_mismatch", cols = c("gray85", "darkred"), pt.size = 0.6) +
  ggtitle("Classification Discrepancies Cells") +
  theme(plot.title = element_text(face = "bold", size = 14))

Idents(query) <- "rf_model_Labels"
p2 <- DimPlot(query, reduction = "umap", label = TRUE, repel = TRUE) + 
  ggtitle("Supervised RF Labels") + 
  theme(legend.position = "none")


# side-by-side in color
options(repr.plot.width = 12, repr.plot.height = 6)
final_plot = p1 + p2 
ggsave( "/Users/kiyounghan/Desktop/Discrepancies_Random_Forest.png", plot = final_plot,width  = 12,height = 10,dpi = 300)

print("==============================================================================")
print(" 13. Downstream Diagnostics...COMPLETED                                       ")
print("==============================================================================")


