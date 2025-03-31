#!/usr/bin/env Rscript

rm(list=ls())
args = commandArgs(trailingOnly=TRUE)
print(args)
message("Number of input arguments: ",  length(args))

if (length(args)!=3){
	stop("Incorrect input argiments. sample name, path to the seurat object containing seeker data, and path to the RDS object containing single cell reference need to be supplied as arguments", call.=FALSE)
}else{
	sample_id = args[1]
	path_to_seurat_object = args[2]
	path_to_reference = args[3]
}

library(spacexr)
library(Matrix)
library(Seurat)
library(ggplot2)

##* FUNCTIONS
# generate output files for Seurat_default or RCTD
process_Seeker_data <- function(sample_name, seurat_object, reference, UMI_min = 0, UMI_min_sigma = 300, gene_count_min = 0, max_cores = 8) {
	message("Importing data")
	sobj <- readRDS(seurat_object) #load seurat object
	count_matrix <- (sobj[["RNA"]]@counts) #get matrix using seurat object
	nUMI <- colSums(count_matrix) #get nUMI - list of total counts
	coords <- as.data.frame(sobj[["SPATIAL"]]@cell.embeddings) #get spatial_coordinates using seurat object
	spatial_RNA <- SpatialRNA(coords, count_matrix, nUMI)
	
	message("Kicking off RCTD")
	RCTD <- create.RCTD(spatial_RNA, readRDS(path_to_reference), max_cores, CELL_MIN_INSTANCE = 0, UMI_min = UMI_min, UMI_min_sigma = UMI_min_sigma)
	RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")
	
	message("Organizing RCTD output")
	RCTD_clusters <- as.data.frame(RCTD@results$results_df["first_type"])
	confidence_levels <- as.data.frame(RCTD@results$results_df["spot_class"])
	df_annotation <- create_annotate_file(as.data.frame(RCTD@results$results_df))
	df_seurat <- data.frame(RCTD@results$results_df["first_type"], RCTD@results$results_df["second_type"], RCTD@results$results_df["spot_class"])
	sobj  <- AddMetaData(object = sobj , metadata = df_seurat[Cells(sobj),], col.name=colnames(df_seurat))
	plot_spatial<- DimPlot(sobj, reduction = 'SPATIAL', group.by = 'first_type')&NoAxes()&theme(plot.title = element_blank())&coord_fixed()
	return(list(RCTD, df_annotation, sobj, plot_spatial))
}

# create annotation file
create_annotate_file <- function(df) {
	df_anno <- data.frame(rownames(df), df["first_type"], df["second_type"], df["spot_class"])
	colnames(df_anno) <- c("Barcode", "First_Cell_Type", "Second_Cell_Type", "Confidence_Level")
	return(df_anno)
}

message("-----------------------------------------")
message("Starting RCTD analysis")

process_data <- process_Seeker_data (sample_name = sample_id, seurat_object = path_to_seurat_object, reference = path_to_reference, gene_count_min = 0, max_cores = 8)

message("RCTD completed")
message("-----------------------------------------")

message("-----------------------------------------")
message("Generating output")

saveRDS(process_data[[1]], paste0(sample_id, "_", "RCTD", ".rds" ))
write.table(process_data[[2]], paste0(sample_id, "_", "RCTD", "_", "annotation.txt"), sep="\t", row.names = FALSE,col.names = TRUE, quote = FALSE)
saveRDS(process_data[[3]], paste0(sample_id, "_", "RCTD", "_", "seurat.rds"))
png(paste0(sample_id, "_", "RCTD", "_", "spatial.png"), units="in", width=12, height=20, res=1200)
plot(process_data[[4]])
dev.off()

message("Output generation completed")
message("-----------------------------------------")
















