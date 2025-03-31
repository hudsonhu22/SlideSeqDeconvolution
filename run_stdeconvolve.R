library(STdeconvolve)

## load  data
obj <- readRDS("/home/hudsonhu/scratch/B13_seurat.rds")
pos <- obj[["SPATIAL"]]@cell.embeddings
cd <- obj[["RNA"]]

## remove pixels with too few genes
counts <- cleanCounts(cd, min.lib.size = 100)

## feature select for genes
corpus <- restrictCorpus(counts, removeAbove=1.0, removeBelow = 0.05)

## choose optimal number of cell-types
ldas <- fitLDA(t(as.matrix(corpus)), Ks = 10)

## get best model results
optLDA <- optimalModel(models = ldas, opt = "min")

## extract deconvolved cell-type proportions (theta) and transcriptional profiles (beta)
results <- getBetaTheta(optLDA, perc.filt = 0.05, betaScale = 1000)
deconProp <- results$theta
deconGexp <- results$beta