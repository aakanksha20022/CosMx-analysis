library("Seurat")
library("tidyverse")
library("glmGamPoi")
library("harmony")
library("SingleR")
library("celldex")
library("devEMF")
library("ggplot2")

#loading in the data--------------
seurat_obj= readRDS("resegmented_cosmx_gut_with_samples.rds")
View(seurat_obj@meta.data)

#renaming the sample column with actual patient_ids-------
#creating a new column for status
sample_list <- list(
  "1" = c("INC2655", "Fuso_neg", "C1"), "2" = c("INC2842", "Fuso_neg", "C1"), "3" = c("INC0223", "Fuso_neg", "C1"), 
  "4" = c("INC0315", "Fuso_neg", "C1"), "5" = c("INC0510", "Fuso_neg", "C1"), "6" = c("INC1662", "Fuso_neg", "C1"), 
  "7" = c("INC0032", "Fuso_neg", "C1"), "8" = c("INC0046", "Fuso_neg", "C1"), "9" = c("INC0074", "Fuso_neg", "C1"), 
  "10" = c("INC2655", "Fuso_pos", "C2"), "11" = c("INC2842", "Fuso_pos", "C2"), "12" = c("INC0223", "Fuso_pos", "C2"), 
  "13" = c("INC0315", "Fuso_pos", "C2"), "14" = c("INC0510", "Fuso_pos", "C2"), "15" = c("INC1662", "Fuso_pos", "C2"), 
  "16" = c("INC0032", "Fuso_pos", "C2"), "17" = c("INC0046", "Fuso_pos", "C2"), "18" = c("INC0074", "Fuso_pos", "C2"), 
  "19" = c("INC2655", "Fuso_pos", "C3"), "20" = c("INC2842", "Fuso_pos", "C3"), "21" = c("INC0223", "Fuso_pos", "C3"), 
  "22" = c("INC0315", "Fuso_pos", "C3"), "23" = c("INC0510", "Fuso_pos", "C3"), "24" = c("INC1662", "Fuso_pos", "C3"),
  "25" = c("INC0027", "Fuso_pos", "C3"), "26" = c("INC0032", "Fuso_pos", "C3"), "27" = c("INC0046", "Fuso_pos", "C3"), 
  "28" = c("INC0074", "Fuso_pos", "C3")
)

# Ensuring the 'sample' column is in character format
seurat_obj@meta.data$sample <- as.character(seurat_obj@meta.data$sample)

# Ensuring 'status' column exists before assigning values
seurat_obj@meta.data$status <- NA

#Additional column to identify each core
seurat_obj@meta.data$core_id= NA

# Looping through each row and update status values
for (i in 1:nrow(seurat_obj@meta.data)) {
  key <- seurat_obj@meta.data$sample[i]  
  if (key %in% names(sample_list)) {
    seurat_obj@meta.data$status[i] <- sample_list[[key]][2] 
  }
}

# Looping through each row and update core numbers
for (i in 1:nrow(seurat_obj@meta.data)) {
  key <- seurat_obj@meta.data$sample[i]  
  if (key %in% names(sample_list)) {
    seurat_obj@meta.data$core_id[i] <- sample_list[[key]][3] 
  }
}

# Looping through each row and update sample ids
for (i in 1:nrow(seurat_obj@meta.data)) {
  key <- seurat_obj@meta.data$sample[i]  
  if (key %in% names(sample_list)) {
    seurat_obj@meta.data$sample[i] <- sample_list[[key]][1]  
  }
}

#additional column to identify pos/neg samples
seurat_obj@meta.data$sample_id <- paste0(seurat_obj@meta.data$sample, "_", seurat_obj@meta.data$status)

#additional column to identify unique samples
seurat_obj@meta.data$unique_id <- paste0(seurat_obj@meta.data$sample, "_", seurat_obj@meta.data$core_id)

#Normalize and scale data
options(future.globals.maxSize = 3 * 1024^3)# Sets limit to 3GB
seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
seurat_obj= ScaleData(seurat_obj)

#calculate variable features
seurat_obj= FindVariableFeatures(seurat_obj)

#PCA
seurat_obj = RunPCA(object = seurat_obj, assay ="RNA", reduction.name="pca_unintegrated")
ElbowPlot(seurat_obj, ndims = 50, reduction = "pca_unintegrated")
dims_to_use =40

# cluster and UMAP
seurat_obj = FindNeighbors(object = seurat_obj, reduction = "pca_unintegrated", dims = 1:dims_to_use)
seurat_obj = FindClusters(object = seurat_obj, reduction = "pca_unintegrated", resolution = 1)
seurat_obj = RunUMAP(seurat_obj, reduction = "pca_unintegrated", reduction.name = "umap_unintegrated", dims = 1:dims_to_use, min.dist = 0.3, spread = 1)

#Visulazation of clusters
DimPlot(seurat_obj, reduction = "umap_unintegrated", group.by = "seurat_clusters", label = TRUE, repel = TRUE) + ggtitle("UMAP by Cluster")
DimPlot(seurat_obj, reduction = "umap_unintegrated", group.by = "sample_id", label=T, repel=T) + ggtitle("UMAP by Sample")
ImageDimPlot(seurat_obj, fov="FusoB.INC.TMA..Sara.10.11.24") + NoLegend()
ImageFeaturePlot(seurat_obj, fov = "FusoB.INC.TMA..Sara.10.11.24" , features = "nFeature_RNA", axes = TRUE)

# Integrate - using harmony
seurat_obj = RunHarmony(seurat_obj,
                        assay.use = "RNA",  # The assay containing gene expression data
                        group.by.vars = c("sample"),  # Metadata column that represents batch effects 
                        theta = c(2),  # Strength of batch correction; typically set between 1-5 (higher values mean stronger correction)
                        reduction = "pca_unintegrated",  # The PCA reduction to use as input for Harmony
                        reduction.save = "pca_integrated")  # Name of the new Harmony-corrected PCA reduction

# re-cluster after integration
seurat_obj = FindNeighbors(seurat_obj, reduction = "pca_integrated", dims = 1:dims_to_use)
seurat_obj = FindClusters(seurat_obj, reduction = "pca_integrated", resolution = 0.5)
seurat_obj = RunUMAP(seurat_obj, dims = 1:dims_to_use, reduction = "pca_integrated", reduction.name ="umap_integrated")

#saving the integrated object
saveRDS(seurat_obj, "seurat_obj_integrated.rds")

seurat_obj= readRDS("seurat_obj_integrated.rds")
DimPlot(seurat_obj, reduction = "umap_integrated", group.by = "seurat_clusters", label = TRUE, repel = TRUE) + ggtitle("UMAP by Cluster")
DimPlot(seurat_obj, reduction = "umap_integrated", group.by = "sample_id", label=T, repel=T) + ggtitle("UMAP by Sample")
# Spatial plot for integrated data
ImageDimPlot(seurat_obj, fov = "FusoB.INC.TMA..Sara.10.11.24", group.by = "seurat_clusters") + NoLegend()

# Spatial feature plot for integrated data
ImageFeaturePlot(seurat_obj, fov = "FusoB.INC.TMA..Sara.10.11.24", features = "nFeature_RNA", axes = TRUE)

#Annotation with SingleR
#loading in the reference
ref= HumanPrimaryCellAtlasData()

#Extracting normalized expressions and converting it to a matrix
seurat_matrix= GetAssayData(seurat_obj, slot= 'data')

cluster_labels= seurat_obj@meta.data$seurat_clusters

singleR_res= SingleR(test=seurat_matrix,
                     ref=ref,
                     labels= ref$label.main,
                     clusters=cluster_labels)

seurat_obj@meta.data$SingleR_annot=singleR_res$labels[match(cluster_labels, rownames(singleR_res))]
DimPlot(seurat_obj, reduction = "umap_integrated", group.by = "SingleR_annot", label = TRUE, repel = TRUE) + ggtitle("UMAP with SingleR Annotations")

#Identifying the top markers
markers= FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, 
                        logfc.threshold = 0.25, test.use = "wilcox") 

top.markers= markers %>%
  group_by(cluster) %>%
  slice_max(n = 15, order_by = avg_log2FC)

write.csv(top.markers, file = "top_markers.csv", row.names = TRUE)
write.csv(markers, file = "markers.csv", row.names = TRUE)

#Heatmap
Idents(seurat_obj)= "seurat_clusters"
clusters_to_keep= as.character(0:12)
seurat_subset= subset(seurat_obj, idents= clusters_to_keep)

#maxcells= min(table(Idents(seurat)))
maxcells= 2000
seurat_subset= subset(seurat_subset, downsample= maxcells)

#filtering for top markers in the top 13 clusters
top.markers_subset <- top.markers %>% 
  filter(cluster %in% clusters_to_keep) %>%  # Keep markers only for selected clusters
  group_by(cluster) %>% 
  slice_max(n = 10, order_by = avg_log2FC)  # Use top 10 markers per cluster

png("Markers_heatmap.png", height = 4000, width = 4000, res = 300)
DoHeatmap(seurat_subset, features = top.markers_subset$gene, size = 4) + 
  scale_fill_gradientn(colors = c("blue", "white", "yellow")) +  # Better contrast
  theme_minimal() +  
  theme(
    axis.text.x = element_blank(),  
    axis.ticks.x = element_blank(),  
    axis.text.y = element_text(size = 8),  
    legend.position = "right"
  ) +
  ggtitle("Top Markers Heatmap")
dev.off()


#Image plot with singleR annotations
png("SpatialPlot_SingleR.png", width = 6000, height = 6000, res = 300)
ImageDimPlot(seurat_obj, fov = "FusoB.INC.TMA..Sara.10.11.24", group.by = "SingleR_annot") + 
  theme(legend.position = "right")
dev.off()

# Rename clusters
seurat_obj = RenameIdents(seurat_obj, 
                             '0' = "Intestinal Epithelial", 
                             '1' = "Epithelial", 
                             '2' = "Immune_T_NK",
                             '3' = "Plasma_B",
                             '4' = "Macrophage_Dendritic_Monocytes",
                             '5' = "Epithelial_Inflam",
                             '6' = "smooth_muscle_fibroblasts",
                             '7' = "Myo_fibroblasts & Stromal",
                             '8' = "B",
                             '9' = "All",
                             '10' = "Endothelial",
                             '11' = "Mast",
                             '12' = "Neuron_Glial",
                             '13' = "Epithelial",
                             '14' = "Epithelial",
                             '15' = "Epithelial",
                             '16' = "Epithelial_Inflam",
                             '17' = "Epithelial_Inflam",
                             '18' = "Epithelial_Inflam",
                             '19' = "Epithelial_Inflam",
                             '20' = "Epithelial_Inflam",
                             '21' = "Epithelial_Inflam")

DimPlot(seurat_obj, label = TRUE, repel = T)

#empty list to store pseudobulk matrices for each cell type
pseudobulk_list <- list()

#pseudobulk for each single R cell type
for (cell_type in unique(Idents(seurat_obj))) {
  seurat_pb = subset(seurat_obj, idents = cell_type)
  print(ncol(seurat_pb))
  pseudobulk = AggregateExpression(seurat_pb, group.by = "sample_id", return.seurat = TRUE)
  pseudobulk_df = data.frame(pseudobulk@assays$RNA$counts)
  pseudobulk_list[[cell_type]] = pseudobulk_df
  write.csv(pseudobulk_df, file = paste0(cell_type, "_pseudobulk.csv"), row.names = TRUE)
}

unique(Idents(seurat_obj))

colnames(seurat_obj@meta.data)
SpatialFeaturePlot(seurat_obj, features = "nCount_RNA")
SpatialDimPlot(seurat_obj)
head(rownames(seurat_obj@meta.data))  
head(rownames(seurat_obj@images$FusoB.INC.TMA..Sara.10.11.24@boundaries))  

Idents(seurat_obj)
#DRAWING SPATIAL REGIONS
# get centroids
centroids <- data.frame(
  x = seurat_obj$x_slide_mm,
  y = seurat_obj$y_slide_mm,
  cluster = as.numeric(as.factor(Idents(seurat_obj))),
  sample = seurat_obj$unique_id
)


names(centroids) = c("X","Y", "CELL_TYPE", "SAMPLE")

# save the centroids
write.table(centroids, "/home/ac715y/scripts/centroids1.csv", quote=FALSE, sep = "\t", col.names = TRUE)

View((seurat_obj@meta.data))

unique(seurat_obj@meta.data$unique_id)
