library(Seurat)
library(spacexr)
library(Matrix)

output_path <- "/home/hudsonhu/scratch/"

# First address the reference
    # need counts, cell_types, nUMI
refObj <- readRDS("/home/hudsonhu/scratch/FreshInstall/RatAllCTOuts/202503211_BC9/trekker_BC9/output/BC9_ConfPositioned_seurat_spatial.rds")
counts <- refObj[["SCT"]]@counts
cell_types <- refObj$seurat_clusters
reference <- Reference(counts, cell_types)

#Save the reference object
saveRDS(reference, file.path(output_path, "ref.rds"))