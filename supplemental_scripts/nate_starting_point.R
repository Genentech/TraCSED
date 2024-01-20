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







#######################################################
##### subset metadata to bc of interest

# time point and treatment
barvar <- 'GFPBC_libB_92196'
bc_interest_metadata <- cell_metadata[cell_metadata$gfp_bc == barvar, ]

# plot the pathways for barvar
for (pathway in names(pathway_genes)) {
  cat("\n\n### ", pathway, "\n")
  pathway_plot <- ggplot(bc_interest_metadata, aes(x=time_point, y=bc_interest_metadata[,paste0('Score.', pathway)], fill=time_point)) +
    geom_violin() + theme_classic() + labs(title=paste0(barvar, paste0(' - ', pathway)),x="time_point", y = paste0('Score.', pathway)) #+ geom_boxplot(width=0.1)
  print(pathway_plot)
}

# plot pathways and treatment for barvar
for (pathway in names(pathway_genes)) {
  cat("\n\n### ", pathway, "\n")
  pathway_plot <- ggplot(bc_interest_metadata, aes(x=time_point, y=bc_interest_metadata[,paste0('Score.', pathway)], fill=treatment)) +
    geom_violin() + theme_classic() + labs(title=paste0(barvar, paste0(' - ', pathway)),x="time_point", y = paste0('Score.', pathway)) #+ geom_boxplot(width=0.1)
  print(pathway_plot)
}





######################################
# old attempt
bar_trx_means = data.frame()
for (bar in levels(barcodes_filtered)){
  if (!(grepl(',', bar)) & grepl('GFPBC', bar)){
    bc_interest_metadata <- cell_metadata_single_bc[cell_metadata_single_bc$gfp_bc == bar, ]
    split_trx <- lapply(split(bc_interest_metadata, ~ treatment, drop = TRUE), split, ~ time_point, drop = TRUE)
    if (length(split_trx) < 1)
      break
    for (trx in 1:length(split_trx)){
      if (length(split_trx[[trx]]) >= 4){
        print(bar)
        #print(length(split_trx[[trx]]))
        subtrx <- split_trx[[trx]]
        for (timp in 1:length(subtrx)){
          trx_bc_paths <- subtrx[[timp]][ , 25:length(subtrx[[timp]])]
          df <- as.data.frame(trx_bc_paths)
          bar_trx_means <- rbind(bar_trx_means, colMeans(df))
          row.names(bar_trx_means)[nrow(bar_trx_means)] <- paste0(bar, paste0(' - ' ,paste0(names(split_trx)[trx], ifelse(names(subtrx)[timp] == 'Month 6', 'Day 180', names(subtrx)[timp]))))
        }
      }
    }
  }
}
colnames(bar_trx_means) <- colnames(bc_interest_metadata[ , 25:length(bc_interest_metadata)])

###### scale pathways across time
combo_trx <- bar_trx_means[rownames(bar_trx_means) %like% 'Palbociclib', ]
combo_trx <- combo_trx[rownames(combo_trx) %like% 'GDC-9545', ]
combo_trx <- scale(combo_trx)
Palbociclib_trx <- bar_trx_means[rownames(bar_trx_means) %like% 'Palbociclib', ]
Palbociclib_trx <- Palbociclib_trx[!(rownames(Palbociclib_trx) %like% 'GDC-9545'), ]
Palbociclib_trx <- scale(Palbociclib_trx)
GDC_trx <- bar_trx_means[rownames(bar_trx_means) %like% 'GDC-9545', ]
GDC_trx <- GDC_trx[!(rownames(GDC_trx) %like% 'Palbociclib'), ]
GDC_trx <- scale(GDC_trx)

# heatmap
heatmap(combo_trx, cexRow = 0.45, cexCol = 0.35)
heatmap(Palbociclib_trx, cexRow = 0.45, cexCol = 0.35)
heatmap(GDC_trx, cexRow = 0.45, cexCol = 0.35)
# complex heatmap

barmatMaker <- function(treatmat){
  barcode <- strsplit(row.names(treatmat), ' - ')
  barcode2 <- strsplit(sapply(barcode, "[[", 2), 'Day ')
  barmat <- as.data.frame(matrix(unlist(barcode),ncol=2,byrow=T))
  barmat <- cbind(barmat, as.data.frame(matrix(unlist(barcode2),ncol=2,byrow=T)))
  colnames(barmat) <- c('Barcode', 'Trash', 'Treatment', 'Day')
  barmat
}


b_combo <- barmatMaker(combo_trx)
b_gdc <- barmatMaker(GDC_trx)
b_palbo <- barmatMaker(Palbociclib_trx)
b_all <- barmatMaker(bar_trx_means)
##########################################
# variables for code below
trx_use <- combo_trx
b_use <- b_combo
###############


rownames(trx_use) <- b_use$Barcode
row_day <- ComplexHeatmap::rowAnnotation(Day = b_use$Day, 
              col = list(Day = c("0" = "white", "1" = "red", "4" = "green", 
                                 "8" = "skyblue", "26" = "blue", "180" = "black")))
row_trx <- ComplexHeatmap::rowAnnotation(Treatment = b_use$Treatment)


barResponse <- function(i){
  x <- barcode_survival[[b_use$Barcode[i]]]
  y <- x[x$treatment == b_use$Treatment[i], ]
  z <- as.numeric(y[y$time_point == paste0('Day ', b_use$Day[i]), 6])
  z
}
SurvStat <- c()
for (i in 1:length(b_use$Barcode)){
  print(b_use$Barcode[i])
  SurvStat[i] <- barResponse(i)
}
row_surv <- ComplexHeatmap::rowAnnotation(Survival = ComplexHeatmap::anno_barplot(SurvStat, 
                            gp = grid::gpar(fill = ifelse(SurvStat >= 0, 'red', 'blue'), lwd = 0)))

ComplexHeatmap::Heatmap(trx_use, name = 'Column Scaled Acitvity', 
                        right_annotation = c(row_day, row_trx, row_surv), 
                        row_names_gp = grid::gpar(fontsize = 4), 
                        column_names_gp = grid::gpar(fontsize = 9))


#x <- barcode_survival[[b_combo$Barcode[1]]]
#y <- x[x$treatment == b_combo$Treatment[1], ]
#z <- as.numeric(y[y$time_point == paste0('Day ', b_combo$Day[1]), 6])

###################################################
# work on this next - fix output so time is cols and barcode-path-activity is row
###################################################
bar_trx_means = data.frame()
time_progression_means = data.frame(matrix(ncol=5,nrow=0, dimnames=list(NULL, c("Day.0", "Day.1", "Day.4", "Day.8", "Day.26"))))
for (bar in levels(barcodes_filtered)){
  if (!(grepl(',', bar)) & grepl('GFPBC', bar)){
    bc_interest_metadata <- cell_metadata_single_bc[cell_metadata_single_bc$gfp_bc == bar, ]
    split_trx <- lapply(split(bc_interest_metadata, ~ treatment, drop = TRUE), split, ~ time_point, drop = TRUE)
    if (length(split_trx) < 1)
      break
    for (trx in 1:length(split_trx)){
      if ('Baseline' %in% names(split_trx)){
        subtrx <- split_trx[[1]]
        pathtime0 <- data.frame()
        trx_bc_paths <- subtrx[[1]][ , 25:length(subtrx[[1]])]
        df <- as.data.frame(trx_bc_paths)
        #bar_trx_means <- rbind(bar_trx_means, colMeans(df))
        #row.names(bar_trx_means)[nrow(bar_trx_means)] <- paste0(bar, paste0(' - ' ,paste0(names(split_trx)[trx], names(subtrx)[timp])))
        pathtime0 <- rbind(pathtime0, colMeans(df))
        colnames(pathtime0) <- colnames(trx_bc_paths)
      }
      if (length(split_trx[[trx]]) >= 4 ){
        print(bar)
        subtrx <- split_trx[[trx]]
        pathtime <- data.frame()
        for (timp in 1:length(subtrx)){
          trx_bc_paths <- subtrx[[timp]][ , 25:length(subtrx[[timp]])]
          df <- as.data.frame(trx_bc_paths)
          #bar_trx_means <- rbind(bar_trx_means, colMeans(df))
          #row.names(bar_trx_means)[nrow(bar_trx_means)] <- paste0(bar, paste0(' - ' ,paste0(names(split_trx)[trx], names(subtrx)[timp])))
          # column name to append to - sub(" ", ".", names(subtrx))[timp]
          #ptime_col <- sub(" ", ".", names(subtrx))[timp]
          pathtime <- rbind(pathtime, colMeans(df))
          colnames(pathtime) <- colnames(trx_bc_paths)
        }
        #pathtime2 <- t(pathtime)
        if (dim(pathtime)[1] >= 4){
          pathtime2 <- rbind(pathtime0, pathtime)
          pathtime2 <- t(pathtime2)
          pathtime2 <- pathtime2[ , 1:5]
          colnames(pathtime2) <- c("Day.0", "Day.1", "Day.4", "Day.8", "Day.26")
          row.names(pathtime2) <- paste0(bar, paste0(' - ' ,paste0(names(split_trx)[trx], colnames(trx_bc_paths))))
          time_progression_means <- rbind(time_progression_means, pathtime2) 
        }
      }
    }
  }
}
#colnames(bar_trx_means) <- colnames(bc_interest_metadata[ , 25:length(bc_interest_metadata)])

# separate
combo_trx <- time_progression_means[rownames(time_progression_means) %like% 'Palbociclib', ]
combo_trx <- combo_trx[rownames(combo_trx) %like% 'GDC-9545', ]

Palbociclib_trx <- time_progression_means[rownames(time_progression_means) %like% 'Palbociclib', ]
Palbociclib_trx <- Palbociclib_trx[!(rownames(Palbociclib_trx) %like% 'GDC-9545'), ]

GDC_trx <- time_progression_means[rownames(time_progression_means) %like% 'GDC-9545', ]
GDC_trx <- GDC_trx[!(rownames(GDC_trx) %like% 'Palbociclib'), ]

# plot pathways and treatment for barvar
trx_df <- combo_trx
treat_var <- 'GDC-9545 + Palbociclib'



MakeRowAttr <- function(treatmat){
  barcode <- strsplit(row.names(treatmat), ' - ')
  barmat <- as.data.frame(matrix(unlist(barcode),ncol=2,byrow=T))
  colnames(barmat) <- c('Barcode', 'Trash')
  barmat$Barcode
}
bar26day <- function(i, treat_var){
  x <- barcode_survival[[trx_rows[i]]]
  y <- x[x$treatment == treat_var, ]
  as.numeric(y[y$time_point == 'Day 26', 6])
}

for (pathway in names(pathway_genes)) {
  cat("\n\n### ", pathway, "\n")
  path_trx <- trx_df[rownames(trx_df) %like% pathway, ]
  trx_rows <- MakeRowAttr(path_trx)
  SurvStat <- c()
  for (i in 1:length(trx_rows)){
    SurvStat[i] <- bar26day(i, treat_var)
  }
  row_surv <- ComplexHeatmap::rowAnnotation(Day26_Survival = ComplexHeatmap::anno_barplot(SurvStat, 
                      gp = grid::gpar(fill = ifelse(SurvStat >= 0, 'red', 'blue'), lwd = 0)))
  # normalize to baseline mean
  pathmean <- path_trx / mean(path_trx$Day.0)
  #pathmean <- scale(pathmean)
  hp <- ComplexHeatmap::Heatmap(as.matrix(pathmean), name = pathway, right_annotation = row_surv,
                          row_names_gp = grid::gpar(fontsize = 4), 
                          column_names_gp = grid::gpar(fontsize = 9),
                          cluster_columns = FALSE, row_labels = trx_rows)
  print(hp)
}


###############################################################################
# plot barcode survival
plot_barcode_set <- c('GFPBC_libB_92196', 'GFPBC_libB_90850', 'GFPBC_libB_37798')

for (barcode in plot_barcode_set) {
  cat("\n\n### ", barcode, "\n")
  barcode_plot <- plot_bcode_freq_tpt(barcode)
  print(barcode_plot)
}

###############################################################################
# output data
surv <- cell_metadata_single_bc[cell_metadata_single_bc$gfp_bc %in% names(barcode_survival), ]
surv <- surv[, c(2,8,9,26:length(surv))]
write.csv(surv, 'pathway_activity_scores.csv')
barsurv <- toJSON(barcode_survival)
###############################################################################
#
###############################################################################
# PLSR
set.seed(1)
model <- plsr(slingPseudotime_Avg~., data=bar_trx_means, scale=TRUE, validation='CV')
summary(model)
plot(model, plottype = "scores", comps = 1:3, )



##############################################################################
assayNames(sce)
dim(logcounts(sce))
sce2 <- sce[ ,colData(sce)$gfp_bc %in% names(barcode_survival)]
dim(logcounts(sce2))
## Remove lowly expressed genes which have less than 10 cells with any counts
sce2 <- sce2[rowSums(counts(sce2) > 1) >= 10, ]
dim(logcounts(sce2))
groups <- colData(sce2)[ ,c('gfp_bc')]
exp <- aggregate(t(counts(sce2)), by=groups, FUN=sum)




##############################################################################
# output the survival data and expression for use in the transformer
######### filter the data
dim(sce)
# only cells in our survival list
sce <- sce[ ,colData(sce)$gfp_bc %in% names(barcode_survival)]
dim(sce)
# only cells in our survival list
sce <- sce[ ,colData(sce)$treatment == 'GDC-9545 200nm']
dim(sce)
# filter so that at least 25% of cells have 3 counts of included genes
sce <- sce[rowSums(counts(sce) > 3) > ncol(sce)/4, ]
dim(sce)

sce <- sce[, colData(sce)$gfp_bc == 'GFPBC_libB_92196']
dim(sce)

bar26day <- function(bar, treat_var, time_var){
  x <- barcode_survival[[bar]]
  y <- x[x$treatment == treat_var, ]
  as.numeric(y[y$time_point == time_var, 6])
}

# change the treatment variable
surv <- c()
treatment1 <- 'GDC-9545 200nm'
for (t in sce$time_point){
  surv <- append(surv, bar26day('GFPBC_libB_92196', treatment1, t))
  print(t)
}
length(surv)
colData(sce)$Survival <- surv

df2 <- t(as.data.frame(logcounts(sce)))
df2 <- as.data.frame(df2)
x <- cbind(sce$slingPseudotime_Avg, cbind(sce$time_point, sce$Survival))
df2 <- cbind(x, df2)
colnames(df2)[1:3] <- c('Pseudotime', 'RealTime', 'Survival')
df2 <- df2[df2$RealTime != 1, ]
df2 <- df2[order( df2[,2], df2[,1] ),]
write.csv(df2,file='bar92196_survexp.csv', row.names=FALSE)

