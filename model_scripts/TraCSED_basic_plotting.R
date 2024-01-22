setwd('/gstore/project/hr_brca_heterogeneity/T47D_trace_Seq_v2')
sce <- gp.sa.core::readResult("traceseq_pathway_score-results/sce", ".", rooted=TRUE)

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
library(entropy)
# added by me
library(ggplot2)
library(data.table)
library(jsonlite)
library(pls)
library(slingshot)

BPPARAM <- MulticoreParam(8)
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
top_barcodes_num <- 20
add_barcodes <- NULL

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
sce <- logNormCounts(sce, BPPARAM = BPPARAM)

########################################################
#
########################################################
#### prepare cell data
cell_metadata <- base::as.data.frame(colData(sce))
cell_metadata$gfp_bc <- fct_explicit_na(cell_metadata$gfp_bc, "")
cell_metadata_single_bc <- cell_metadata %>%
  rownames_to_column() %>%
  filter(gfp_bc_count == "1")

##### prepare cell condition data
cell_treatments <- levels(cell_metadata$treatment)
cell_time_points <- levels(cell_metadata$time_point)
baseline_time_point <- cell_time_points[1]
final_time_point <- tail(cell_time_points, n = 1)

##### my additions
cell_barcode_groups <- levels(cell_metadata$gfp_bc)

##### helper functions
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

calculate_entropy <- function(x) {
  return(entropy(table(x)))
}

##### prepare barcode statistics data
barcode_cluster_summary <- cell_metadata_single_bc %>%
  group_by(label) %>%
  mutate(cluster_cell_count = n()) %>%
  ungroup() %>%
  group_by(gfp_bc, label) %>%
  mutate(barcode_cluster_cell_count = n()) %>%
  ungroup() %>%
  mutate(
    barcode_cluster_frequency = barcode_cluster_cell_count / 
      cluster_cell_count
  ) %>%
  group_by(gfp_bc) %>%
  summarise(
    max_cluster_frequency = max(barcode_cluster_frequency) 
  )
entropy_summary <- cell_metadata_single_bc %>%
  group_by(gfp_bc, time_point) %>%
  summarise(entropy = calculate_entropy(label)) %>%
  mutate(
    time_point = paste0("entropy_", gsub(" ", "_", tolower(time_point)))
  ) %>%
  spread(time_point, entropy, fill = NA)
barcode_summary <- cell_metadata_single_bc %>%
  group_by(gfp_bc) %>%
  summarise(
    cluster_count = n_distinct(label), 
    treatment_count = n_distinct(treatment), 
    time_point_count = n_distinct(time_point),
    overall_frequency = n() / nrow(cell_metadata_single_bc), 
  ) %>%
  left_join(barcode_cluster_summary) %>%
  left_join(entropy_summary) %>%
  as.data.frame()

####### function for barcode ploting ########
#
#############################################
all_cell_data <- cell_metadata_single_bc %>%
  group_by(treatment, time_point) %>%
  summarise(all_cell_count = n())

barcode_survival <- list()
for (barcode in barcode_summary$gfp_bc){
  barcode_cell_data <- cell_metadata_single_bc %>%
    filter(gfp_bc == barcode) %>%
    group_by(treatment, time_point) %>%
    summarise(barcode_cell_count = n()) %>%
    full_join(all_cell_data) %>%
    replace_na(list(barcode_cell_count = 0)) %>%
    mutate(barcode_frequency = barcode_cell_count / all_cell_count)
  baseline <- as.numeric(barcode_cell_data[barcode_cell_data$treatment == 'Baseline', 5])
  barcode_cell_data$barcode_response <- log(barcode_cell_data$barcode_frequency / baseline, 2)
  b <- list(barcode_cell_data)
  barcode_survival <- append(barcode_survival, b)
  print(barcode)
}
names(barcode_survival) <- barcode_summary$gfp_bc

plot_bcode_freq_tpt <- function(barcode) {
  all_cell_data <- cell_metadata_single_bc %>%
    group_by(treatment, time_point) %>%
    summarise(all_cell_count = n())
  barcode_cell_data <- cell_metadata_single_bc %>%
    filter(gfp_bc == barcode) %>%
    group_by(treatment, time_point) %>%
    summarise(barcode_cell_count = n()) %>%
    full_join(all_cell_data) %>%
    replace_na(list(barcode_cell_count = 0)) %>%
    mutate(barcode_frequency = barcode_cell_count / all_cell_count)
  ggplot(barcode_cell_data,
         mapping = aes(x = time_point,
                       y = barcode_frequency,
                       fill = treatment)) +
    geom_col(position = "dodge2") +
    labs(title = "Barcode frequency by time point and treatment",
         subtitle = glue("Barcode: {barcode}"),
         x = "Time point",
         y = "Barcode frequency",
         fill = "Treatment") +
    theme_bw()
}

bar26day <- function(bar, treat_var){
  x <- barcode_survival[[bar]]
  y <- x[x$treatment == treat_var, ]
  as.numeric(y[y$time_point == 'Day 26', 6])
}

##### potential barcodes of interest
barcodes_filtered <- barcode_summary %>%
  filter(
    treatment_count == length(levels(cell_metadata$treatment)),
    time_point_count == length(levels(cell_metadata$time_point))
  ) %>%
  slice_max(max_cluster_frequency, n = top_barcodes_num) %>%
  pull(gfp_bc)
barcodes_of_interest <- barcode_summary %>%
  filter(gfp_bc %in% c(barcodes_filtered, add_barcodes))
knitr::kable(barcodes_of_interest[, c(1, 2, 5, 6)])
knitr::kable(barcodes_of_interest[, c(1, 7:ncol(barcodes_of_interest))])


#### Dimensionality reduction
# run umap
set.seed(seed_value)
sce <- runUMAP(sce, dimred = "PCA", BPPARAM = BPPARAM)

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
# Cluster
dittoDimPlot(sce, "label", reduction.use = "UMAP", size = point_size)


# individual day treatment
for (t in levels(sce$time_point)) {
  cat("\n\n### ", t, "\n")
  umap_plot <- dittoDimPlot(sce[ ,colData(sce)$time_point == t], "treatment", 
                            reduction.use = "UMAP", size = point_size, main = t)
  print(umap_plot)
}
# run umap for individual days
set.seed(seed_value)
for (t in levels(sce$time_point)){
  cat("\n\n### ", t, "\n")
  iso <- c(t, "Day 0")
  sct <- runUMAP(sce[ ,colData(sce)$time_point %in% iso], dimred = "PCA", BPPARAM = BPPARAM)
  sct <- sct[ ,colData(sct)$treatment != 'GDC-9545 + Palbociclib']
  umap_plot <- dittoDimPlot(sct[ ,colData(sct)$time_point %in% iso], "treatment", 
                            reduction.use = "UMAP", size = point_size, main = paste0(t, ' - Time Isolated UMAP Run'))
  print(umap_plot)
  #pca_plot <- dittoDimPlot(sct[ ,colData(sct)$time_point %in% iso], "treatment", 
  #                         reduction.use = "PCA", size = point_size, main = paste0(t, ' - Time Isolated PCA Run'))
  #print(pca_plot)
}

keepbar <- c('GFPBC_libB_92196', 'GFPBC_libB_20182')
sce_red <- sce[ ,colData(sce)$gfp_bc %in% keepbar]
sce_red <- sce_red[ ,colData(sce_red)$time_point != 'Month 6']
sce_palbo <- sce_red[ ,colData(sce_red)$treatment == 'Palbociclib 200nm']
sce_9545 <- sce_red[ ,colData(sce_red)$treatment == 'GDC-9545 200nm']
sce_base <- sce_red[ ,colData(sce_red)$treatment == 'Baseline']
umap_plot <- dittoDimPlot(sce_palbo, "gfp_bc", reduction.use = "UMAP", size = point_size, 
                          split.by = "time_point", legend.show = TRUE, main = 'Palbociclib')
print(umap_plot)

umap_plot <- dittoDimPlot(sce_9545, "gfp_bc", reduction.use = "UMAP", size = point_size, 
                          split.by = "time_point", legend.show = TRUE, main = 'GDC-9545')
print(umap_plot)

umap_plot <- dittoDimPlot(sce_base, "gfp_bc", reduction.use = "UMAP", size = point_size, 
                          split.by = "time_point", legend.show = TRUE, main = 'Baseline')
print(umap_plot)


day_9545 <- sce_base[ ,colData(sce_base)$time_point == 'Day 0']
day_9545_92196 <- day_9545[ ,colData(day_9545)$gfp_bc == 'GFPBC_libB_92196']
day_9545_20182 <- day_9545[ ,colData(day_9545)$gfp_bc == 'GFPBC_libB_20182']
mean(reducedDim(day_9545_92196, 'UMAP')[, 1]) # mean x
mean(reducedDim(day_9545_92196, 'UMAP')[, 2]) # mean y
mean(reducedDim(day_9545_20182, 'UMAP')[, 1]) # mean x
mean(reducedDim(day_9545_20182, 'UMAP')[, 2]) # mean y
#u <- reducedDim(day1_9545_92196, 'UMAP')

#######################
# attempt to visualize barcodes
# takes some time
umap_plot <- dittoDimPlot(sce, "gfp_bc", reduction.use = "UMAP", size = point_size, 
                          split.by = "time_point", legend.show = FALSE)
print(umap_plot)

#######################################################################
# cluster barcodes
#######################################################################
single_bcs <- unique(sce$gfp_bc)
single_bcs <- as.data.frame(single_bcs[!(grepl(',', single_bcs))])
colnames(single_bcs) <- "single_barcodes"
sbc <- sce[ ,colData(sce)$gfp_bc %in% single_bcs$single_barcodes]

for (t in levels(sce$time_point)) {
  cat("\n\n### ", t, "\n")
  pca_plot <- dittoDimPlot(sbc[ ,colData(sbc)$time_point == t], "gfp_bc", 
                           reduction.use = "PCA", size = point_size, main = t, 
                           legend.show = FALSE)
  #png(filename=paste0(t, "pca.png"), width=900,)
  par(las=3)
  pca_data <- pca_plot$plot_env$Target_data
  group_ordered <- with(pca_data, reorder(color, X, median))
  boxplot(X ~ group_ordered, pca_data, cex.axis = 0.15, xlab = NULL, 
          ylab = 'PC1', main = t)
  #plot(pca_data$color, pca_data$X, cex.axis = 0.5, xlab = NULL, ylab = 'PC1')
  print(pca_plot)
  #dev.off()
}


# umap for individual days
for (t in levels(sce$time_point)) {
  cat("\n\n### ", t, "\n")
  umap_plot <- dittoDimPlot(sbc[ ,colData(sbc)$time_point == t], labels, 
                            reduction.use = "UMAP", size = point_size, main = t)
  #legend.show = FALSE)
  #png(filename=paste0(t, "pca.png"), width=1200, height=550)
  par(las=3)
  umap_data <- umap_plot$plot_env$Target_data
  group_ordered <- with(umap_data, reorder(color, X, median))
  boxplot(X ~ group_ordered, umap_data, cex.axis = 0.15, xlab = NULL, 
          ylab = 'UMAP1', main = t)
  print(umap_plot)
  #dev.off()
}



# UMAP of just GFPBC_libB_92196
sbc <- sce[ ,colData(sce)$gfp_bc == 'GFPBC_libB_92196']
umap_plot <- dittoDimPlot(sbc, "treatment", reduction.use = "UMAP", size = point_size, 
                          split.by = "time_point", legend.show = TRUE, main = 'GFPBC_libB_92196')
print(umap_plot)


# parent umap for individual days
# labels can be gfp_bc or labels param below
t <- 'Day 0'
sbc2 <- sbc[ ,colData(sbc)$time_point == t]
sbc2 <- sbc2[ ,colData(sbc2)$gfp_bc %in% barcode_summary$gfp_bc]
labels <- c()
for (bar in sbc2$gfp_bc){
  labels <- append(labels, ifelse(bar26day(bar, "GDC-9545 + Palbociclib") >= 0, 'Growth', 'Depleted'))
}
umap_plot <- dittoDimPlot(sbc2, labels, reduction.use = "UMAP", 
                          size = point_size, main = t)
print(umap_plot)

###################################
# umap for individual days
# labels can be gfp_bc or labels param below
set.seed(seed_value)
t <- 'Day 0'
sbc2 <- sbc[ ,colData(sbc)$time_point == t]
sbc2 <- sbc2[ ,colData(sbc2)$gfp_bc %in% barcode_summary$gfp_bc]
sbc2 <- runUMAP(sbc2, dimred = "PCA", BPPARAM = BPPARAM)
labels <- c()
for (bar in sbc2$gfp_bc){
  labels <- append(labels, ifelse(bar26day(bar, "GDC-9545 + Palbociclib") >= 0, 'Growth', 'Depleted'))
}
umap_plot <- dittoDimPlot(sbc2, labels, reduction.use = "UMAP", 
                          size = point_size, main = t)
print(umap_plot)




#######################################################################
#
#######################################################################
# pathway params
pathway_file <- "~/git/HTSEQ-1022-analysis-pipeline-for-trace-seq/data/pathway_genes_NGS4327.tab"
gene_id_column <- "symbol"
point_size <- 0.1

# prep geneset data
pathway_df <- read.delim(pathway_file, header = FALSE)
pathway_genes <- strsplit(pathway_df[, 2], ",")
names(pathway_genes) <- pathway_df[, 1]

# calculate pathway score
#pathway_scores <- sumCountsAcrossFeatures(sce, pathway_genes,exprs_values = "logcounts",average = TRUE)
#pathway_scores <- as.data.frame(t(pathway_scores))
#colnames(pathway_scores) <- paste0("Score.", colnames(pathway_scores))

# add pathway scores to the SCE object
#colData(sce) <- cbind(colData(sce), pathway_scores)

#### violin plots
# time point and treatment
for (pathway in names(pathway_genes)) {
  cat("\n\n### ", pathway, "\n")
  pathway_plot <- dittoPlot(sce[,colData(sce)$treatment != 'GDC-9545 + Palbociclib'], paste0("Score.", pathway),
                            group.by = "time_point",
                            color.by = "treatment",
                            plots = "vlnplot", vlnplot.lineweight = 0.5)
  print(pathway_plot)
}

#### ER signature
er_sig <- sce$Score.ER_induced - sce$Score.ER_repressed
pathway_plot <- dittoPlot(sce, er_sig,
                          group.by = "time_point",
                          color.by = "treatment",
                          plots = "vlnplot", vlnplot.lineweight = 0.5)
print(pathway_plot)
pathway_plot <- dittoPlot(sce, er_sig, group.by = "label", plots = "vlnplot", 
                          vlnplot.lineweight = 0.5, ylab = 'induced - repressed',
                          main = 'ER signature')
print(pathway_plot)

# cluster
for (pathway in names(pathway_genes)) {
  cat("\n\n### ", pathway, "\n")
  pathway_plot <- dittoPlot(sce, paste0("Score.", pathway),
                            group.by = "label", legend.show = FALSE,
                            plots = "vlnplot", vlnplot.lineweight = 0.5)
  print(pathway_plot)
}

# cluster and treatment
for (pathway in names(pathway_genes)) {
  cat("\n\n### ", pathway, "\n")
  pathway_plot <- dittoPlot(sce, paste0("Score.", pathway),
                            group.by = "label", color.by = "treatment",
                            plots = "vlnplot", vlnplot.lineweight = 0.5)
  print(pathway_plot)
}
