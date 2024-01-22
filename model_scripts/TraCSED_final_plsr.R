setwd('/gstore/project/hr_brca_heterogeneity/T47D_trace_Seq_v2')
sce <- gp.sa.core::readResult("traceseq_pathway_score-results/sce", ".", rooted=TRUE)

# Palbociclib is resolution0.75
# GDC-9545 is resolution1
# Combination is resolution1.25


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
library(plspm)
library(plsdepot)
library(DESeq2)
library(biomaRt)
library(org.Hs.eg.db)
library(fgsea)
library(genomitory)
library(Seurat)
library(batchelor)
library(proxy)
library(ggridges)

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

# pathway params
pathway_file <- "~/git/HTSEQ-1022-analysis-pipeline-for-trace-seq/data/pathway_genes_NGS4327.tab"
gene_id_column <- "symbol"
point_size <- 0.1

# prep geneset data
pathway_df <- read.delim(pathway_file, header = FALSE)
pathway_genes <- strsplit(pathway_df[, 2], ",")
names(pathway_genes) <- pathway_df[, 1]

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
# summary stats
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
  barcode_cell_data$barcode_response <- barcode_cell_data$barcode_frequency / baseline
  b <- list(barcode_cell_data)
  barcode_survival <- append(barcode_survival, b)
  print(barcode)
}
names(barcode_survival) <- barcode_summary$gfp_bc

bar26day <- function(bar, treat_var, time_var){
  x <- barcode_survival[[bar]]
  y <- x[x$treatment == treat_var, ]
  as.numeric(y[y$time_point == time_var, 6])
}


######################################

# PLSR
set.seed(1)
######### filter the data
dim(sce)
# only cells in our survival list
sce <- sce[ ,colData(sce)$gfp_bc %in% names(barcode_survival)]
dim(sce)
# only cells at day 0
sce <- sce[ ,colData(sce)$time_point == 'Day 0']
dim(sce)
# use this filter because stringency is necessary here but not for the model
# filter so that at least 25% of cells have 3 counts of included genes
sce <- sce[rowSums(counts(sce) > 3) > ncol(sce)/4, ]
dim(sce)

#################################################################
######### get the survival at time_point
bar26day <- function(bar, treat_var, time_var){
  x <- barcode_survival[[bar]]
  y <- x[x$treatment == treat_var, ]
  as.numeric(y[y$time_point == time_var, 6])
}

# change the treatment variable
day26surv1 <- c()
day8surv1 <- c()
day4surv1 <- c()
day1surv1 <- c()
day26surv2 <- c()
day8surv2 <- c()
day4surv2 <- c()
day1surv2 <- c()
day26surv3 <- c()
day8surv3 <- c()
day4surv3 <- c()
day1surv3 <- c()
treatment1 <- 'Palbociclib 200nm'
treatment2 <- 'GDC-9545 200nm'
treatment3 <- 'GDC-9545 + Palbociclib'
for (bar in sce$gfp_bc){
  day26surv1 <- append(day26surv1, bar26day(bar, treatment1, 'Day 26'))
  day8surv1 <- append(day8surv1, bar26day(bar, treatment1, 'Day 8'))
  day4surv1 <- append(day4surv1, bar26day(bar, treatment1, 'Day 4'))
  day1surv1 <- append(day1surv1, bar26day(bar, treatment1, 'Day 1'))
  day26surv2 <- append(day26surv2, bar26day(bar, treatment2, 'Day 26'))
  day8surv2 <- append(day8surv2, bar26day(bar, treatment2, 'Day 8'))
  day4surv2 <- append(day4surv2, bar26day(bar, treatment2, 'Day 4'))
  day1surv2 <- append(day1surv2, bar26day(bar, treatment2, 'Day 1'))
  day26surv3 <- append(day26surv3, bar26day(bar, treatment3, 'Day 26'))
  day8surv3 <- append(day8surv3, bar26day(bar, treatment3, 'Day 8'))
  day4surv3 <- append(day4surv3, bar26day(bar, treatment3, 'Day 4'))
  day1surv3 <- append(day1surv3, bar26day(bar, treatment3, 'Day 1'))
  print(bar)
}
length(day26surv3)
colData(sce)$Day26Survival <- day26surv3

# model data
df2 <- t(as.data.frame(logcounts(sce)))
df2 <- as.data.frame(df2)
x <- day26surv3
model <- plsr(x ~ ., data=df2, ncomp=30, scale=TRUE, validation='CV')
summary(model)
plot(model, plottype = "scores", comps = 1:5)


###### pca ######
set.seed(seed_value)
sce <- runPCA(sce, BPPARAM = BPPARAM)
sce <- runUMAP(sce, dimred = "PCA", BPPARAM = BPPARAM)
colData(sce)$UMAP_1 <- reducedDim(sce, 'UMAP')[, 1]
colData(sce)$UMAP_2 <- reducedDim(sce, 'UMAP')[, 2]

# ridge distances between resistant barcodes
posbars <- c('GFPBC_libB_92196', 'GFPBC_libB_37798')
scepos <- sce[,colData(sce)$gfp_bc %in% posbars]
df <- as.data.frame(colData(scepos)[, c("UMAP_1", "UMAP_2")])
posdistmat <- dist(df, df, method = 'euclidean')
posdist <- as.data.frame(as.numeric(posdistmat))
lab <- rep("Positive/Positive", dim(posdist)[1])
posdist <- cbind(posdist, lab)
colnames(posdist) <- c('Distance', 'Label')

# ridge distances between sensitive barcodes
negbars <- c('GFPBC_libB_20182', 'GFPBC_libB_31581', 'GFPBC_libB_40758', 
             'GFPBC_libB_94397', 'GFPBC_libB_83026', 'GFPBC_libB_81397', 
             'GFPBC_libB_10678', 'GFPBC_libB_90850', 'GFPBC_libB_10737')
sceneg <- sce[,colData(sce)$gfp_bc %in% negbars]
df2 <- as.data.frame(colData(sceneg)[, c("UMAP_1", "UMAP_2")])
negdistmat <- dist(df2, df2, method = 'euclidean')
negdist <- as.data.frame(as.numeric(negdistmat))
lab <- rep("Negative/Negative", dim(negdist)[1])
negdist <- cbind(negdist, lab)
colnames(negdist) <- c('Distance', 'Label')

# ridge distances from resistant to sensitive
alldistmat <- dist(df, df2, method = 'euclidean')
alldist <- as.data.frame(as.numeric(alldistmat))
lab <- rep("Positive/Negative", dim(alldist)[1])
alldist <- cbind(alldist, lab)
colnames(alldist) <- c('Distance', 'Label')

totaldf <- rbind(posdist, rbind(negdist, alldist))

# ridge plot for PCA
ggplot(totaldf, aes(x = Distance, y = Label, fill = Label)) +
  geom_density_ridges() +
  theme_ridges() + 
  theme(legend.position = "none")

################## plsr
dim(reducedDim(sce, "PCA"))
plsr_scores <- scores(model)[,1:5] ### change slice to get more components
reducedDim(sce, "PLSR") <- plsr_scores
dim(reducedDim(sce, "PLSR"))

# run umap
set.seed(1)
sce <- runUMAP(sce, dimred = "PLSR", BPPARAM = BPPARAM)
colData(sce)$UMAP_1 <- reducedDim(sce, 'UMAP')[, 1]
colData(sce)$UMAP_2 <- reducedDim(sce, 'UMAP')[, 2]

# ridge distances between resistant barcodes
posbars <- c('GFPBC_libB_92196', 'GFPBC_libB_37798')
scepos <- sce[,colData(sce)$gfp_bc %in% posbars]
df <- as.data.frame(colData(scepos)[, c("UMAP_1", "UMAP_2")])
posdistmat <- dist(df, df, method = 'euclidean')
posdist <- as.data.frame(as.numeric(posdistmat))
lab <- rep("Positive/Positive", dim(posdist)[1])
posdist <- cbind(posdist, lab)
colnames(posdist) <- c('Distance', 'Label')

# ridge distances between sensitive barcodes
negbars <- c('GFPBC_libB_20182', 'GFPBC_libB_31581', 'GFPBC_libB_40758', 
             'GFPBC_libB_94397', 'GFPBC_libB_83026', 'GFPBC_libB_81397', 
             'GFPBC_libB_10678', 'GFPBC_libB_90850', 'GFPBC_libB_10737')
sceneg <- sce[,colData(sce)$gfp_bc %in% negbars]
df2 <- as.data.frame(colData(sceneg)[, c("UMAP_1", "UMAP_2")])
negdistmat <- dist(df2, df2, method = 'euclidean')
negdist <- as.data.frame(as.numeric(negdistmat))
lab <- rep("Negative/Negative", dim(negdist)[1])
negdist <- cbind(negdist, lab)
colnames(negdist) <- c('Distance', 'Label')

# ridge distances from resistant to sensitive
alldistmat <- dist(df, df2, method = 'euclidean')
alldist <- as.data.frame(as.numeric(alldistmat))
lab <- rep("Positive/Negative", dim(alldist)[1])
alldist <- cbind(alldist, lab)
colnames(alldist) <- c('Distance', 'Label')

totaldf <- rbind(posdist, rbind(negdist, alldist))

# ridge plot for PLSR
ggplot(totaldf, aes(x = Distance, y = Label, fill = Label)) +
  geom_density_ridges() +
  theme_ridges() + 
  theme(legend.position = "none")

####################################################
# clustering
cluster_resolutions <- c(0.1, 0.2, 0.3, 0.5, 0.75, 1, 1.25, 1.5, 2)

set.seed(seed_value)
snn_graph <- buildSNNGraph(sce, use.dimred = "PLSR", k = k_value,
                           BPPARAM = BPPARAM, BNPARAM = BNPARAM)
cluster_modularity_list <- list()
for (cluster_resolution in cluster_resolutions){
  set.seed(seed_value)
  cell_clusters <- cluster_leiden(
    snn_graph, resolution_parameter = cluster_resolution
  )$membership
  sce[[glue("resolution{cluster_resolution}")]] <- as.factor(cell_clusters)
  cluster_modularity_list[[glue("resolution{cluster_resolution}")]] <-
    pairwiseModularity(snn_graph, cell_clusters, as.ratio = TRUE)
}

clustree(sce, prefix = "resolution") +
  guides(edge_alpha = FALSE, edge_colour = FALSE)

for (cluster_resolution in cluster_resolutions) {
  cat("\n\n### ", glue("resolution={cluster_resolution}"), "\n")
  umap_plot <- dittoDimPlot(sce, glue("resolution{cluster_resolution}"),
                            reduction.use = "UMAP",
                            size = point_size, do.label = TRUE,
                            labels.highlight = FALSE, legend.show = TRUE)
  print(umap_plot)
}

scater::plotUMAP(sce, colour_by = 'resolution1.25')

# gene names
ids <- row.names(sce)
annots <- select(org.Hs.eg.db, keys=ids, 
                 columns="SYMBOL", keytype="ENSEMBL")
annots <- annots[!duplicated(annots$ENSEMBL), ]
row.names(annots) <- annots$ENSEMBL
annots$GENE <- ifelse(is.na(annots$SYMBOL), annots$ENSEMBL, annots$SYMBOL)
rownames(sce) <- annots$GENE

# markers
cluster.markers.any <- findMarkers(sce, groups = sce$resolution1.25, 
                                   pval.type = "any")
cluster.markers.all <- findMarkers(sce, groups = sce$resolution1.25, 
                                   pval.type = "all")

# top 40 gene heatmap
marker_genes <- unique(rownames(cluster.markers.all$`3`))[1:40]
heatmap_plot <- plotGroupedHeatmap(sce, marker_genes, group = "resolution1.25", 
                                   center = TRUE, cex=0.85)

##### add in the correct resolutions
heat <- scater::plotHeatmap(sce, features = unique(unlist(lapply(cluster.markers.any, function(w) rownames(w)[1:5]))),
                            columns = colnames(sce)[order(sce$resolution1.25)],   ##### add in the correct resolutions
                            colour_columns_by = "resolution1.25", cluster_cols = FALSE,
                            show_colnames = FALSE, cluster_rows = FALSE)
genes <- unique(rownames(cluster.markers.all$`3`))[1:10]
for (g in genes){
  print(g)
  v <- dittoPlot(sce, g, plots = "vlnplot", group.by = "resolution1.25", vlnplot.lineweight = 0.5)
  print(v)
  u <- scater::plotUMAP(sce, colour_by = g)
  print(u)
}

genes <- c('SNHG25', 'PLK2', 'COX1', 'COX2', 'CLDN1', 'MGP', 'CRIP1', 'CRIP2', 'SNCG')
for (g in genes){
  print(g)
  v <- dittoPlot(sce, g, plots = "vlnplot", group.by = "resolution1.25", vlnplot.lineweight = 0.5)
  print(v)
  u <- scater::plotUMAP(sce, colour_by = g)
  print(u)
}

##########################################################
# heatmap of cluster assignments and barcode abundance
assignment <- as.data.frame(cbind(as.data.frame(sce$gfp_bc), sce$resolution1.25))
colnames(assignment) <- c("barcode", "cluster")
df <- as.data.frame(table(assignment))
df2 <- reshape(df, timevar = "cluster", idvar = "barcode", direction = "wide")
row.names(df2) <- df2[,1]
df2 <- df2[,-1]
# remove any barcode with less than 5 cells
df2 <- df2[rowSums(df2[])>5,]
barcode_heatmap <- df2/rowSums(df2)
#barcode_heatmap <- t(scale(t(df2)))


bar26day <- function(i, treat_var, time_var){
  x <- barcode_survival[[row.names(barcode_heatmap)[i]]]
  y <- x[x$treatment == treat_var, ]
  log(as.numeric(y[y$time_point == time_var, 6]), 2)
}
SurvStat1 <- c()
SurvStat2 <- c()
SurvStat3 <- c()
SurvStat4 <- c()
treat <- 'GDC-9545 + Palbociclib'
for (i in 1:length(row.names(df2))){
  print(row.names(barcode_heatmap)[i])
  sur <- bar26day(i, treat, 'Day 1')
  SurvStat1[i] <- ifelse(sur != -Inf, sur, -3)
  sur <- bar26day(i, treat, 'Day 4')
  SurvStat2[i] <- ifelse(sur != -Inf, sur, -3)
  sur <- bar26day(i, treat, 'Day 8')
  SurvStat3[i] <- ifelse(sur != -Inf, sur, -3)
  sur <- bar26day(i, treat, 'Day 26')
  SurvStat4[i] <- ifelse(sur != -Inf, sur, -3)
}
row_surv1 <- ComplexHeatmap::rowAnnotation(D1 = ComplexHeatmap::anno_barplot(SurvStat1, 
                                                                             gp = grid::gpar(fill = ifelse(SurvStat1 >= 0, 'red', 'blue'), lwd = 0)))
row_surv2 <- ComplexHeatmap::rowAnnotation(D4 = ComplexHeatmap::anno_barplot(SurvStat2, 
                                                                             gp = grid::gpar(fill = ifelse(SurvStat2 >= 0, 'red', 'blue'), lwd = 0)))
row_surv3 <- ComplexHeatmap::rowAnnotation(D8 = ComplexHeatmap::anno_barplot(SurvStat3,
                                                                             gp = grid::gpar(fill = ifelse(SurvStat3 >= 0, 'red', 'blue'), lwd = 0)))
row_surv4 <- ComplexHeatmap::rowAnnotation(D26 = ComplexHeatmap::anno_barplot(SurvStat4,
                                                                              gp = grid::gpar(fill = ifelse(SurvStat4 >= 0, 'red', 'blue'), lwd = 0)))


ComplexHeatmap::Heatmap(barcode_heatmap, name = 'Fraction', 
                        right_annotation = c(row_surv1, row_surv2, row_surv3, row_surv4), 
                        row_names_gp = grid::gpar(fontsize = 4), 
                        column_names_gp = grid::gpar(fontsize = 9), 
                        column_labels = c("C1", "C2", "C3", "C4", "C5", "C6"))

##########################################################################