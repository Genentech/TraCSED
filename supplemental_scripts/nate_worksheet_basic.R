setwd('/gstore/project/hr_brca_heterogeneity/T47D_trace_Seq_v2')
#sce <- gp.sa.core::readResult("traceseq_pathway_score-results/sce", ".", rooted=TRUE)
sce <- gp.sa.core::readResult("traceseq_qc-results/sce", rooted=TRUE)

library(tidyverse)
library(glue)
library(scales)
library(viridis)
library(knitr)
library(formattable)
library(gp.sa.core)
library(BiocParallel)
library(BiocNeighbors)
library(Matrix)
library(DelayedArray)
library(SingleCellExperiment)
library(scater)
library(scran)
library(PCAtools)
library(igraph)
library(bluster)
library(dittoSeq)
library(clustree)
library(pheatmap)

BPPARAM <- MulticoreParam(4)
BNPARAM <- AnnoyParam()

#### analysis params
assay_to_use <- "counts_all"
cell_min_genes <- 500
cell_max_percent_mt <- 25
species <- "human"
cc_markers <- readRDS(system.file("exdata",glue("{species}_cycle_markers.rds"), package = "scran"))
top_hvgs_num <- 2000
k_value <- 50
cluster_resolutions <- c(0.05, 0.1, 0.15, 0.2)
point_size <- 0.1
seed_value <- 42

#### choose the assay for further analysis
counts(sce) <- as(assay(sce, assay_to_use), "dgCMatrix")

#### cell filtering
# filter cells by number of expressed genes
sce <- sce[, sce$detected >= cell_min_genes]
ncol(sce)
#filter cells by percentage of mitochondrial reads
sce <- sce[, sce$subsets_mito_percent <= cell_max_percent_mt]
ncol(sce)
# stats
summary(sce$sum)
summary(sce$detected)
table(sce$sample)
table(sce$library)
table(sce$treatment)
table(sce$time_point)

#### Data normalization
########## These commented lines don't work ################
quick_clusters <- quickCluster(sce, block = sce$library, block.BPPARAM = BPPARAM) 
#sce <- computeSumFactors(sce, clusters = quick_clusters, BPPARAM = BPPARAM)
sce <- logNormCounts(sce, BPPARAM = BPPARAM)

#### cell cycle annotation
set.seed(seed_value)
cyclone_results <- cyclone(sce, cc_markers, BPPARAM = BPPARAM)
sce$phase <- as.factor(cyclone_results$phases)
# summary
table(sce$phase)
# detect cell cycle-related genes to be removed
diff <- getVarianceExplained(sce, "phase")
is_cc_related <- diff > 5
summary(is_cc_related)

#### feature selection
# model gene variance
dec <- modelGeneVar(sce, BPPARAM = BPPARAM)
# plot gene variance
dec %>%
  as.data.frame() %>%
  ggplot(mapping = aes(x = mean, y = total)) +
  geom_point() +
  geom_smooth() +
  labs(
    title = "Modelling gene variance",
    x = "Mean log-expression",
    y = "Total variance"
  ) +
  theme_bw()
# select top highly variable genes (HVGs)
top_hvgs <- getTopHVGs(dec[which(!is_cc_related),], n = top_hvgs_num)
length(top_hvgs)

#### Dimensionality reduction
# run the PCA and remove PCs corresponding to noise
set.seed(seed_value)
sce <- denoisePCA(sce, technical = dec, subset.rowm = top_hvgs, BPPARAM = SerialParam())
#number of significant PCs
ncol(reducedDim(sce, "PCA"))
# scree plot elbow
percent_var <- attr(reducedDim(sce, "PCA"), "percentVar")
pca_elbow <- findElbowPoint(percent_var)
pca_elbow
# scree plot
ggplot(mapping = aes(x = seq_along(percent_var), y = percent_var)) +
  geom_point() +
  geom_vline(aes(xintercept = ncol(reducedDim(sce, "PCA")),
                 colour = "Denoised PCs")) +
  geom_vline(aes(xintercept = pca_elbow,
                 colour = "Elbow point")) +
  scale_colour_manual(values = c("Denoised PCs" = "red",
                                 "Elbow point" = "blue")) +
  labs(
    title = "Scree plot",
    x = "PC",
    y = "Variance explained (%)"
  ) +
  theme_bw()

# run umap
set.seed(seed_value)
sce <- runUMAP(sce, dimred = "PCA") # , BPPARAM = BPPARAM)

#### UMAP plots
# Sample
dittoDimPlot(sce, "sample", reduction.use = "UMAP", size = point_size)
# Library
dittoDimPlot(sce, "library", reduction.use = "UMAP", size = point_size)
# Treatment
dittoDimPlot(sce, "treatment", reduction.use = "UMAP", size = point_size)
dittoDimPlot(sce, "treatment", reduction.use = "UMAP", size = point_size,
             split.by = "time_point", split.ncol = 2)
# time point
dittoDimPlot(sce, "time_point", reduction.use = "UMAP", size = point_size)
dittoDimPlot(sce, "time_point", reduction.use = "UMAP", size = point_size,
             split.by = "treatment", split.ncol = 2)
# Phase
dittoDimPlot(sce, "phase", reduction.use = "UMAP", size = point_size)
# UMI count
dittoDimPlot(sce, "sum", reduction.use = "UMAP", size = point_size,
             min = quantile(sce$sum, 0.05),
             max = quantile(sce$sum, 0.95))
# Detected gene count
dittoDimPlot(sce, "detected", reduction.use = "UMAP", size = point_size,
             min = quantile(sce$detected, 0.05),
             max = quantile(sce$detected, 0.95))
# MT percent
dittoDimPlot(sce, "subsets_mito_percent", reduction.use = "UMAP",
             size = point_size,
             min = quantile(sce$subsets_mito_percent, 0.05),
             max = quantile(sce$subsets_mito_percent, 0.95))
# GFP barcode count
dittoDimPlot(sce, "gfp_bc_count", reduction.use = "UMAP", size = point_size)
