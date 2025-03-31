library(Seurat)
library(spacexr)
library(Matrix)

# First address the reference
    # need counts, cell_types, nUMI
refObj <- readRDS("/home/hudsonhu/scratch/FreshInstall/RatAllCTOuts/202503211_BC9/trekker_BC9/output/BC9_ConfPositioned_seurat_spatial.rds")
counts <- refObj[["RNA"]]
cell_types <- refObj$seurat_clusters
reference <- Reference(counts, cell_types)

# Now the data to be deconvolved
    # need coords (data.frame), counts (untransformed matrix, row is gene), nUMI (optional)

obj <- readRDS("/home/hudsonhu/scratch/B13_seurat.rds")
coords <- obj[["SPATIAL"]]@cell.embeddings
counts <- obj[["RNA"]]
puck <- SpatialRNA(coords, counts)

# Start running rctd
# myRCTD <- create.RCTD(puck, reference, max_cores = 16)
# results <- myRCTD@results
# norm_weights = normalize_weights(results$weights) # to make them sum to one
# cell_type_names <- myRCTD@cell_type_info$info[[2]] # list of cell type names
# spatialRNA <- myRCTD@spatialRNA
