rm(list=ls())
library(anndata)
library(basilisk)
library(data.table)
library(dplyr)
library(ggplot2)
library(reticulate)
library(Seurat)
library(splatter)
library(zellkonverter)

# Paths are relative to this script's directory (citeseq/R/), which is what
# ZILN.Rproj sets as the R working directory. The arm's committed input lives in
# citeseq/data/; everything this chain produces goes to citeseq/results/.
data_dir <- "../data/"
results_dir <- "../results/"
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)


use_python("/Users/sjun6/opt/anaconda3/envs/nbsr/bin/python")

methods <- c("wilcox", "wilcox_limma", "t", "MAST", "negbinom")
output_path <- results_dir
if (!dir.exists(output_path)) {
  dir.create(output_path, recursive = T)
}

MIN_READS <- 3
MIN_CELLS <- 5

dat <- readRDS(paste0(results_dir, "pbmc10k_cd4_memory.rds"))

cts <- LayerData(dat, assay = "RNA", layer="counts")
dim(cts)
row_idxs <- which(rowSums(cts > MIN_READS) >= MIN_CELLS)
length(row_idxs)

sub_dat <- dat[row_idxs,]
out_path <- paste0(data_dir, "memory_CD4.h5ad")
sce <- as.SingleCellExperiment(sub_dat)
writeH5AD(sce, file = out_path)

