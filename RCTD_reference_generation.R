# R Script for Single-Cell Reference Generation for RCTD

# Load necessary libraries
library(spacexr)

# Set the path to directory with count matrix and meta data files
files_dir <- ""

# Load in counts matrix
# Make sure row.names are gene names
counts <- read.csv(file.path(files_dir, "counts.csv"), row.names=1) 

#Load in meta_data (annotation)
metadata <- read.csv(file.path(files_dir, "metadata.csv")) 
#head(metadata$cells)
#head( metadata$cell.types)

#Create cell_types named list
cell_types <- metadata$cell.types; names(cell_types) <- metadata$cells

#Convert to factor data type
cell_types <- as.factor(cell_types)
#cell_types

#Create the Reference object
RCTD_reference <-  Reference(counts, cell_types)

#Save the reference object
saveRDS(RCTD_reference, file.path(files_dir, "ref.rds"))
