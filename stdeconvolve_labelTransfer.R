library(Seurat)
library(harmony)

# load data
seurat1 <- readRDS("")
seurat2 <- readRDS("")

# 1. Merge Seurat objects
seurat_merged <- merge(seurat1, y = seurat2, add.cell.ids = c("Sample1", "Sample2"))

# 2. Normalize the data
seurat_merged <- NormalizeData(seurat_merged)

# 3. Find variable features
seurat_merged <- FindVariableFeatures(seurat_merged, selection.method = "vst", nfeatures = 2000)

# 4. Scale the data
seurat_merged <- ScaleData(seurat_merged)

# 5. Run PCA
seurat_merged <- RunPCA(seurat_merged, npcs = 30)

# 6. Integrate using Harmony
seurat_merged <- RunHarmony(seurat_merged, group.by.vars = "batch")  # Replace "batch" with your batch column in metadata

# Map the unlabelled cluster to reference labels
cluster_mapping <- c(
  "0" = "T cells",
  "1" = "B cells",
  "2" = "Monocytes",
  "3" = "NK cells"
)

# Plot clusters
seurat_merged$cell_type <- factor(cluster_mapping[as.character(seurat_merged$seurat_clusters)])
DimPlot(seurat_merged, reduction = "umap", group.by = "cell_type")
