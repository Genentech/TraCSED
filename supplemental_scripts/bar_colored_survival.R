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
library(biomaRt)
library(org.Hs.eg.db)
library(ggridges)
library(fgsea)
library(genomitory)

BPPARAM <- MulticoreParam(16)
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




##################################################################
# gene names
# think if adding pathways downstream before adding
ids <- row.names(sce)
annots <- select(org.Hs.eg.db, keys=ids, 
                 columns="SYMBOL", keytype="ENSEMBL")
annots <- annots[!duplicated(annots$ENSEMBL), ]
row.names(annots) <- annots$ENSEMBL
annots$GENE <- ifelse(is.na(annots$SYMBOL), annots$ENSEMBL, annots$SYMBOL)
rownames(sce) <- annots$GENE



##############################################################
######### filter the data
dim(sce)
# only cells in our survival list
sce <- sce[ ,colData(sce)$gfp_bc %in% names(barcode_survival)]
dim(sce)
# only cells in our survival list
treats <- c('GDC-9545 200nm', 'Baseline')
sce <- sce[ ,colData(sce)$treatment %in% treats]
dim(sce)
# only cells up to day 26
times <- c('Day 0', 'Day 1', 'Day 4', 'Day 8', 'Day 26')
sce <- sce[ ,colData(sce)$time_point %in% times]
dim(sce)
# filter so that at least 25% of cells have 3 counts of included genes
sce <- sce[rowSums(counts(sce) > 3) > ncol(sce)/10, ]
#sce <- sce[rowSums(counts(sce) > 1) > ncol(sce)/10, ]
dim(sce)

###############################
# cgp pathways filtered 
cgp_pathways <- getFeatureSetCollection("GMTY42:human/C2.CGP.gmt.bz2@REVISION-1")
pathways_genes2 <- as.list(cgp_pathways)
paths <- as.data.frame(elementMetadata(cgp_pathways)[1])
names(pathways_genes2) <- paths$name
p <- pathways_genes2[grepl('BREAST',names(pathways_genes2))]
p <- append(p, pathways_genes2[grepl('ESR',names(pathways_genes2))])
p <- lapply(p, function(genes) {
  genes <- which(rowData(sce)[['gene_id']] %in% genes)
  if (length(genes) == 0) {
    genes <- NULL
  }
  return(genes)
})
pathway_scores <- sumCountsAcrossFeatures(sce, p,
                                          exprs_values = "logcounts",
                                          average = TRUE)
pathway_scores <- as.data.frame(t(pathway_scores))
colnames(pathway_scores) <- paste0("Score.", colnames(pathway_scores))
colData(sce) <- cbind(colData(sce), pathway_scores)



# pca umap
set.seed(1)
sce <- runUMAP(sce, dimred = "PCA", BPPARAM = BPPARAM)
dittoDimPlot(sce, "time_point", reduction.use = "UMAP", size = point_size)
barsofi <- c('GFPBC_libB_92196', 'GFPBC_libB_90850', 'GFPBC_libB_10678',
             'GFPBC_libB_20182', 'GFPBC_libB_37798', 'GFPBC_libB_31581',
             'GFPBC_libB_40758', 'GFPBC_libB_94397', 'GFPBC_libB_83026')
dittoDimPlot(sce[ ,colData(sce)$gfp_bc %in% barsofi], "gfp_bc", reduction.use = "UMAP", size = 0.3)

dittoDimPlot(sce, "PIP", reduction.use = "UMAP", size = point_size)
dittoDimPlot(sce, "Score.STEIN_ESRRA_TARGETS_UP", reduction.use = "UMAP", size = point_size)

var <- "Score.STEIN_ESRRA_TARGETS_UP"
for (bar in barsofi){
  pathway_plot <- dittoPlot(sce[ ,colData(sce)$gfp_bc == bar], var, group.by = "time_point", plots = "vlnplot", 
                            vlnplot.lineweight = 0.5, ylab = 'ESRRA targets score',
                            main = bar)
  print(pathway_plot)
}


umap_plot <- dittoDimPlot(sce[ ,colData(sce)$gfp_bc == 'GFPBC_libB_90850'], "PIP", reduction.use = "UMAP", size = point_size, 
                          split.by = "time_point", legend.show = TRUE)
print(umap_plot)

umap_plot <- dittoDimPlot(sce[ ,colData(sce)$gfp_bc == 'GFPBC_libB_90850'], "Score.STEIN_ESRRA_TARGETS_UP", reduction.use = "UMAP", size = point_size, 
                          split.by = "time_point", legend.show = TRUE)
print(umap_plot)

# pca trajectory
# clustering
cluster_resolutions <- c(0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.75, 0.85, 1, 1.25)
snn_graph <- buildSNNGraph(sce, use.dimred = "PCA", k = k_value,
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

clustree(sce, prefix = "resolution") + guides(edge_alpha = FALSE, edge_colour = FALSE)

for (cluster_resolution in cluster_resolutions) {
  cat("\n\n### ", glue("resolution={cluster_resolution}"), "\n")
  umap_plot <- dittoDimPlot(sce, glue("resolution{cluster_resolution}"),
                            reduction.use = "UMAP",
                            size = point_size, do.label = TRUE,
                            labels.highlight = FALSE, legend.show = TRUE)
  print(umap_plot)
}

dittoDimPlot(sce, "time_point", reduction.use = "UMAP", size = point_size)

# run slingshot
start_cluster <- "2"
end_cluster <- "8"
sce <- slingshot(sce, reducedDim = "PCA", stretch = 0, approx_points = 200)
sce$slingPseudotime_Global <- sce$slingPseudotime_1
sce$slingPseudotime_1 <- NULL
sce$slingClusters <- NULL
slingPseudotime(sce)

sce <- slingshot(sce, reducedDim = "PCA",
                 clusterLabels = sce$resolution0.75,
                 start.clus = start_cluster,
                 end.clus = end_cluster,
                 stretch = 0, approx_points = 200)
slingPseudotime(sce)
sce$slingPseudotime_Avg2 <- rowMeans(slingPseudotime(sce), na.rm = TRUE)

lineages <- getLineages(data = reducedDim(sce, type = "UMAP"), 
                        clusterLabels = as.numeric(sce$resolution0.75), 
                        start.clus = start_cluster, end.clus = end_cluster) #define where to start the trajectories
branch_assignments <- slingBranchID(sce)
curves <- getCurves(lineages, approx_points = 300, thresh = 0.01, stretch = 0.8, allow.breaks = FALSE, shrink = 0.99)
cellWeights <- slingCurveWeights(curves)
#View(cellWeights)

# Avg
dittoDimPlot(sce, "slingPseudotime_Avg2", reduction.use = "UMAP",
             size = point_size) + scale_color_viridis_c(direction = -1)
# trajectory
dittoDimPlot(sce, "resolution0.3", reduction.use = "UMAP", size = point_size,
             do.label = TRUE, legend.show = FALSE,
             add.trajectory.lineages = slingLineages(sce),
             trajectory.cluster.meta = "resolution0.3")


####### ridge in terms of pnorm
### ridge plot
# as.character vs not is different view
scetime <- as.character(sce$time_point)  #scetime <- as.character(sce$time_point)
scetime <- cbind(scetime, sce$slingPseudotime_Avg2)  #scetime <- cbind(scetime, as.numeric(sce$slingPseudotime_Avg2))
scetime <- cbind(scetime, as.character(sce$gfp_bc))
colnames(scetime) <- c('time_point', 'pseudotime', 'barcode')
scetime <- as.data.frame(scetime)
scetime$pseudotime <- as.numeric(scetime$pseudotime)

dittoPlot(sce, "slingPseudotime_Avg2", group.by = "time_point",
          plots = c("ridgeplot"), ridgeplot.lineweight = 0.5,
          legend.show = FALSE)

###### normalize based on 1% and 99%
###### if below then 0, if above then 1
normalize <- function(x, globalmin, globalmax, na.rm = TRUE) {
  norm <- c()
  for (i in x){
    if (i >= globalmax){
      norm <- append(norm, 1)
    }
    if (i <= globalmin){
      norm <- append(norm, 0)
    }
    if(i < globalmax & i > globalmin){
      norm <- append(norm, (i- globalmin) /(globalmax-globalmin)) 
    } 
  }
  return(norm)
}

# ks test
#gmin <- min(sce$slingPseudotime_Avg2)
#gmax <- max(sce$slingPseudotime_Avg2)
g1 <- quantile(sce$slingPseudotime_Avg2, probs = 0.01)
g99 <- quantile(sce$slingPseudotime_Avg2, probs = 0.99)
pnorm <- normalize(sce$slingPseudotime_Avg2, g1, g99)
sce$Normalized_pseudotime <- pnorm
p0 <- sce[,colData(sce)$time_point == 'Day 0']$Normalized_pseudotime
p1 <- sce[,colData(sce)$time_point == 'Day 1']$Normalized_pseudotime
p4 <- sce[,colData(sce)$time_point == 'Day 4']$Normalized_pseudotime
p8 <- sce[,colData(sce)$time_point == 'Day 8']$Normalized_pseudotime
p26 <- sce[,colData(sce)$time_point == 'Day 26']$Normalized_pseudotime

k10 <- ks.test(p1, p0)
k41 <- ks.test(p4, p1)
k84 <- ks.test(p8, p4)
k268 <- ks.test(p26, p8)
#k268$statistic
k268$statistic + k84$statistic + k41$statistic + k10$statistic


k40 <- ks.test(p4, p0)
k80 <- ks.test(p8, p0)
k260 <- ks.test(p26, p0)
k81 <- ks.test(p8, p1)
k261 <- ks.test(p26, p1)
k264 <- ks.test(p26, p4)
k268$statistic + k84$statistic + k41$statistic + k10$statistic + 
  k40$statistic + k80$statistic + k260$statistic + k81$statistic +
  k261$statistic + k264$statistic


### ridge plot
dittoPlot(sce, "Normalized_pseudotime", group.by = "time_point",
          plots = c("ridgeplot"), ridgeplot.lineweight = 0.5,
          legend.show = FALSE)

# normtime
dittoDimPlot(sce, "Normalized_pseudotime", reduction.use = "UMAP",
             size = point_size) + scale_color_viridis_c(direction = -1)

######## survival cutoff using variance
######## plot variance in survival based on ptime, no cells, all together, and individual

###############
# PLSR
# filter for plsr model
# PLSR
set.seed(1)
# model data
df2 <- t(as.data.frame(logcounts(sce)))
df2 <- as.data.frame(df2)
x <- as.vector(colData(sce)$time_point)
x <- replace(x, x == "Day 0", 1)
x <- replace(x, x == 'Day 1', 2)
x <- replace(x, x == 'Day 4', 3)
x <- replace(x, x == 'Day 8', 4)
x <- replace(x, x == 'Day 26', 5)

# model it
model <- plsr(x ~ ., data=df2, ncomp=30, scale=TRUE, validation='CV')
summary(model)
plot(model, plottype = "scores", comps = 1:5)

# add plsr components
dim(reducedDim(sce, "PCA"))
plsr_scores <- scores(model)[,1:5] ### change slice to get more components
reducedDim(sce, "PLSR") <- plsr_scores
dim(reducedDim(sce, "PLSR"))

###### plsr trajectory ########
set.seed(1)
sce <- runUMAP(sce, dimred = "PLSR", BPPARAM = BPPARAM)

# treatment umap
dittoDimPlot(sce, "treatment", reduction.use = "UMAP", size = point_size)
dittoDimPlot(sce, "time_point", reduction.use = "UMAP", size = point_size)


# clustering
cluster_resolutions <- c(0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.75, 0.8, 0.85, 0.9, 0.95, 1, 1.25)
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

clustree(sce, prefix = "resolution") + guides(edge_alpha = FALSE, edge_colour = FALSE)

for (cluster_resolution in cluster_resolutions) {
  cat("\n\n### ", glue("resolution={cluster_resolution}"), "\n")
  umap_plot <- dittoDimPlot(sce, glue("resolution{cluster_resolution}"),
                            reduction.use = "UMAP",
                            size = point_size, do.label = TRUE,
                            labels.highlight = FALSE, legend.show = TRUE)
  print(umap_plot)
}

# run slingshot
start_cluster <- "1"
end_cluster <- "7"
sce <- slingshot(sce, reducedDim = "PLSR", stretch = 0, approx_points = 200)
sce$slingPseudotime_Global <- sce$slingPseudotime_1
sce$slingPseudotime_1 <- NULL
sce$slingClusters <- NULL
slingPseudotime(sce)

sce <- slingshot(sce, reducedDim = "PLSR",
                 clusterLabels = sce$resolution0.3,
                 start.clus = start_cluster,
                 end.clus = end_cluster,
                 stretch = 0, approx_points = 200)
slingPseudotime(sce)
sce$slingPseudotime_Avg2 <- rowMeans(slingPseudotime(sce), na.rm = TRUE)

lineages <- getLineages(data = reducedDim(sce, type = "UMAP"), 
                        clusterLabels = as.numeric(sce$resolution0.3), 
                        start.clus = start_cluster, end.clus = end_cluster) #define where to start the trajectories
branch_assignments <- slingBranchID(sce)
curves <- getCurves(lineages, approx_points = 300, thresh = 0.01, stretch = 0.8, allow.breaks = FALSE, shrink = 0.99)
cellWeights <- slingCurveWeights(curves)
#View(cellWeights)

# Avg
dittoDimPlot(sce, "slingPseudotime_Avg2", reduction.use = "UMAP",
             size = point_size) + scale_color_viridis_c(direction = -1)
# trajectory
dittoDimPlot(sce, "resolution0.75", reduction.use = "UMAP", size = point_size,
             do.label = TRUE, legend.show = FALSE,
             add.trajectory.lineages = slingLineages(sce),
             trajectory.cluster.meta = "resolution0.75")

### ridge plot
# as.character vs not is different view
scetime <- as.character(sce$time_point)  #scetime <- as.character(sce$time_point)
scetime <- cbind(scetime, sce$slingPseudotime_Avg2)  #scetime <- cbind(scetime, as.numeric(sce$slingPseudotime_Avg2))
scetime <- cbind(scetime, as.character(sce$gfp_bc))
colnames(scetime) <- c('time_point', 'pseudotime', 'barcode')
scetime <- as.data.frame(scetime)
scetime$pseudotime <- as.numeric(scetime$pseudotime)

dittoPlot(sce, "slingPseudotime_Avg2", group.by = "time_point",
          plots = c("ridgeplot"), ridgeplot.lineweight = 0.5,
          legend.show = FALSE)

###### normalize based on 1% and 99%
###### if below then 0, if above then 1
normalize <- function(x, globalmin, globalmax, na.rm = TRUE) {
  norm <- c()
  for (i in x){
    if (i >= globalmax){
      norm <- append(norm, 1)
    }
    if (i <= globalmin){
      norm <- append(norm, 0)
    }
    if(i < globalmax & i > globalmin){
      norm <- append(norm, (i- globalmin) /(globalmax-globalmin)) 
    } 
  }
  return(norm)
}

# ks test
#gmin <- min(sce$slingPseudotime_Avg2)
#gmax <- max(sce$slingPseudotime_Avg2)
g1 <- quantile(sce$slingPseudotime_Avg2, probs = 0.01)
g99 <- quantile(sce$slingPseudotime_Avg2, probs = 0.99)
pnorm <- normalize(sce$slingPseudotime_Avg2, g1, g99)
sce$Normalized_pseudotime <- pnorm
p0 <- sce[,colData(sce)$time_point == 'Day 0']$Normalized_pseudotime
p1 <- sce[,colData(sce)$time_point == 'Day 1']$Normalized_pseudotime
p4 <- sce[,colData(sce)$time_point == 'Day 4']$Normalized_pseudotime
p8 <- sce[,colData(sce)$time_point == 'Day 8']$Normalized_pseudotime
p26 <- sce[,colData(sce)$time_point == 'Day 26']$Normalized_pseudotime

k10 <- ks.test(p1, p0)
k41 <- ks.test(p4, p1)
k84 <- ks.test(p8, p4)
k268 <- ks.test(p26, p8)
#k268$statistic
k268$statistic + k84$statistic + k41$statistic + k10$statistic


k40 <- ks.test(p4, p0)
k80 <- ks.test(p8, p0)
k260 <- ks.test(p26, p0)
k81 <- ks.test(p8, p1)
k261 <- ks.test(p26, p1)
k264 <- ks.test(p26, p4)
k268$statistic + k84$statistic + k41$statistic + k10$statistic + 
  k40$statistic + k80$statistic + k260$statistic + k81$statistic +
  k261$statistic + k264$statistic


### ridge plot
dittoPlot(sce, "Normalized_pseudotime", group.by = "time_point",
          plots = c("ridgeplot"), ridgeplot.lineweight = 0.5,
          legend.show = FALSE)

# normtime
dittoDimPlot(sce, "Normalized_pseudotime", reduction.use = "UMAP",
             size = point_size) + scale_color_viridis_c(direction = -1)




################################################################
# order numerically not by character!!!!
scetime <- scetime[order(scetime[,2]), ]
barsurv <- c()
for (i in 1:length(scetime$barcode)){
  s <- length(which(scetime[c(i:dim(scetime)[1]),]$barcode == scetime$barcode[i])) # number of survivng clones
  stot <- s / length(scetime[c(i:dim(scetime)[1]),]$barcode) # no surv clone / no of all surviving
  d <- length(which(scetime[c(1:i),]$barcode == scetime$barcode[i])) # no of depleted clones
  dtot <- d / length(scetime[c(1:i),]$barcode) # no depl clone / no of all depleted
  print('Surv')
  print(length(scetime[c(i:dim(scetime)[1]),]$barcode)) # no of all surviving
  print(s) # number of survivng clones
  print(stot) # no surv clone / no of all surviving
  print('Depl')
  print(length(scetime[c(1:i),]$barcode)) # no of all depleted
  print(d) # no of depleted clones
  print(dtot) # no depl clone / no of all depl
  alltot <- stot / dtot
  print('All')
  print(alltot)
  barsurv <- append(barsurv, alltot)
}

scetime$barsurvival <- barsurv
##################################################
# plot survival vs time order for barcode
bar <- 'GFPBC_libB_92196' # plsr res0.25 day0 - 0.159 --- day26 - 0.225 - Surv
bar2 <- 'GFPBC_libB_20182' # plsr res0.25 day0 - 0.124 --- day26 - 0.07
bar3 <- 'GFPBC_libB_90850' # plsr res0.25 day0 - 0.080 --- day26 - 0.35 - Surv
bar4 <- 'GFPBC_libB_10678' # plsr res0.25 day0 - 0.044 --- day26 - 0.153 - Surv
bar5 <- 'GFPBC_libB_37798' # plsr res0.25 day0 - 0.129  --- day26 - 0.026
bar6 <- 'GFPBC_libB_31581' # plsr res0.25 day0 - 0.044 --- day26 - 0.020
bar7 <- 'GFPBC_libB_40758' # plsr res0.25 day0 - 0.040  --- day26 - 0.0073
bar8 <- 'GFPBC_libB_94397' # plsr res0.25 day0 - 0.042 --- day26 - 0.013
bar9 <- 'GFPBC_libB_83026' # plsr res0.25 day0 - 0.052 --- day26 - 0.015
bar10 <- 'GFPBC_libB_10737'
bar11 <- 'GFPBC_libB_81397'
scetime0 <- scetime[scetime$time_point == 'Day 0', ]
scetime26 <- scetime[scetime$time_point == 'Day 26', ]

dim(scetime0[scetime0$barcode == bar, ])[1]/dim(scetime0)[1]
dim(scetime26[scetime26$barcode == bar, ])[1]/dim(scetime26)[1]

minMax <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

scetime_bar <- scetime[scetime$barcode == bar2, ]
plot(1:dim(scetime_bar)[1], scetime_bar$barsurvival, main=bar, xlab='Order', ylab='Survival', ylim=c(0, 2.5))


scetime_bar <- scetime[scetime$barcode == bar9, ]
changerate <- c()
for (i in 1:length(scetime_bar$barsurvival)-1){
  changerate <- append(changerate, abs(scetime_bar$barsurvival[i+1] - scetime_bar$barsurvival[i]))
}

plot(scetime_bar$pseudotime[-1], minMax(changerate), type='l', col='red', main='Change over time', xlab='Pseudotime', ylab='Normalized Abs Change', ylim=c(0, 1))
lines(scetime_bar$pseudotime[-1], minMax(changerate), col='yellow')


colors <- c('red', 'blue', 'green', 'black', 'purple', 'yellow', 'pink', 'orange', 'cyan', 'darkred', 'darkblue')
plot((1:length(changerate))/length(changerate), minMax(changerate), type='l', col='red', main='Change over cells', xlab='Fraction of progress', ylab='Normalized Abs Change', ylim=c(0, 1))
lines((1:length(changerate))/length(changerate), minMax(changerate), col='darkblue')

# find the biggest gap between normalized absolute change >= ?
# keep all cells in the biggest gap, remove all cells outside boundary

cellbounds <- function(change, scebar){
  norm <- minMax(change)
  position <- which(norm >= 0.15)
  half <- length(change)/2
  poslow <- max(position[position <= half])
  poshigh <- min(position[position > half])
  scebar <- scebar[poslow:poshigh, ]
  return(scebar)
}

bars <- c(bar, bar2, bar3, bar4, bar5, bar6, bar7, bar8, bar9, bar10, bar11)
for (b in bars){
  scetime_bar <- scetime[scetime$barcode == bar11, ]
  changerate <- c()
  for (i in 1:length(scetime_bar$barsurvival)-1){
    changerate <- append(changerate, abs(scetime_bar$barsurvival[i+1] - scetime_bar$barsurvival[i]))
  }
  scebar <- cellbounds(changerate, scetime_bar)
}



##################################################
# match order of sce object according to barsurv
scetime2 <- scetime[match(rownames(colData(sce)), rownames(scetime)),]
#scetime_bar <- scetime[scetime$barcode == bar7, ]
#plot(1:dim(scetime_bar)[1], scetime_bar$barsurvival, main=bar7, xlab='Order', ylab='Survival', ylim=c(0, 2.5))


sce$barsurvival <- scetime2$barsurvival
# if using a cutoff on day 0
d0_med <- median(sce[ ,colData(sce)$time_point == 'Day 0']$slingPseudotime_Avg2)
d26_q <- quantile(sce[ ,colData(sce)$time_point == 'Day 26']$slingPseudotime_Avg2, probs= c(0.75))
dim(sce)
sce2 <- sce[ ,colData(sce)$slingPseudotime_Avg2 >= d0_med]
dim(sce2)
#ptime_300 <- scetime$pseudotime[length(scetime$pseudotime)-300]
sce2 <- sce2[ ,colData(sce2)$slingPseudotime_Avg2 <= 40] #### was d26_q
dim(sce2)
sce <- sce2

dittoPlot(sce, "slingPseudotime_Avg2", group.by = "time_point",
          plots = c("ridgeplot"), ridgeplot.lineweight = 0.5,
          legend.show = FALSE)



# as.character vs not is different view
scetime3 <- as.character(sce$time_point)  #scetime <- as.character(sce$time_point)
scetime3 <- cbind(scetime3, sce$slingPseudotime_Avg2)  #scetime <- cbind(scetime, as.numeric(sce$slingPseudotime_Avg2))
scetime3 <- cbind(scetime3, as.character(sce$gfp_bc))
colnames(scetime3) <- c('time_point', 'pseudotime', 'barcode')
scetime3 <- as.data.frame(scetime3)
scetime3$pseudotime <- as.numeric(scetime3$pseudotime)


scetime4 <- scetime3[match(rownames(scetime), rownames(scetime3)),]
scetime4$barsurvival <- scetime$barsurvival
scetime4 <- scetime4[complete.cases(scetime4), ]

scetime_bar <- scetime4[scetime4$barcode == bar9, ]
plot(1:dim(scetime_bar)[1], scetime_bar$barsurvival, main=bar9, xlab='Order', ylab='Survival')

bar <- 'GFPBC_libB_92196' # plsr res0.25 day0 - 0.159 --- day26 - 0.225 - Surv
bar2 <- 'GFPBC_libB_20182' # plsr res0.25 day0 - 0.124 --- day26 - 0.07
bar3 <- 'GFPBC_libB_90850' # plsr res0.25 day0 - 0.080 --- day26 - 0.35 - Surv
bar4 <- 'GFPBC_libB_10678' # plsr res0.25 day0 - 0.044 --- day26 - 0.153 - Surv
bar5 <- 'GFPBC_libB_37798' # plsr res0.25 day0 - 0.129  --- day26 - 0.026
bar6 <- 'GFPBC_libB_31581' # plsr res0.25 day0 - 0.044 --- day26 - 0.020
bar7 <- 'GFPBC_libB_40758' # plsr res0.25 day0 - 0.040  --- day26 - 0.0073
bar8 <- 'GFPBC_libB_94397' # plsr res0.25 day0 - 0.042 --- day26 - 0.013
bar9 <- 'GFPBC_libB_83026' #

##################################################
# FOR GENES
# output json with pseudotime surv {Barcode_1: df, Barcode_n: df_n} structure
barleft <- unique(sce$gfp_bc)
barlist <- intersect(barleft, names(barcode_survival))
treats <- c('GDC-9545 200nm', 'Baseline')
tdata <- list()
tnames <- c()
##### turn the order by to ptime only or realtime + ptime
for (bar in barlist){
  bar_sce <- sce[, colData(sce)$gfp_bc == bar]
  if (dim(bar_sce)[2] > 25){   #### this line is the minimum number of cells 
    print(bar)
    print(dim(bar_sce)[2])
    df <- t(as.data.frame(logcounts(bar_sce)))
    df <- as.data.frame(df)
    x <- cbind(bar_sce$slingPseudotime_Avg2, cbind(bar_sce$time_point, bar_sce$barsurvival))
    df <- cbind(x, df)
    colnames(df)[1:3] <- c('Pseudotime', 'RealTime', 'Survival')
    #df <- df[order( df[,2], df[,1] ),] # order by realtime and then order internal by ptime
    df <- df[order(df[,1]), ] # order by pseudotime
    print(dim(df))
    b <- list(df)
    tdata <- append(tdata, b)
    tnames <- append(tnames, bar)
  }
}
names(tdata) <- tnames
tjson <- toJSON(tdata)
write(tjson, "/gstore/project/hr_brca_heterogeneity/T47D_trace_Seq_v2/Nathan/tran_9545_barcodes75_pcares0.25.json")


##################################################
##################################################
# FOR PATHWAYS
# output json with pseudotime surv {Barcode_1: df, Barcode_n: df_n} structure
barleft <- unique(sce$gfp_bc)
barlist <- intersect(barleft, names(barcode_survival))
treats <- c('GDC-9545 200nm', 'Baseline')
tdata <- list()
tnames <- c()
##### turn the order by to ptime only or realtime + ptime
for (bar in barlist){
  bar_sce <- sce[, colData(sce)$gfp_bc == bar]
  if (dim(bar_sce)[2] > 25){   #### this line is the minimum number of cells 
    print(bar)
    print(dim(bar_sce)[2])
    df <- colData(bar_sce)[24:length(colData(bar_sce))]
    df <- df[ ,grepl("Score",colnames(df))]
    x <- cbind(bar_sce$slingPseudotime_Avg2, cbind(bar_sce$time_point, bar_sce$barsurvival))
    df <- cbind(x, df)
    colnames(df)[1:3] <- c('Pseudotime', 'RealTime', 'Survival')
    #df <- df[order( df[,2], df[,1] ),] # order by realtime and then order internal by ptime
    df <- df[order(df[,1]), ] # order by pseudotime
    print(dim(df))
    b <- list(df)
    tdata <- append(tdata, b)
    tnames <- append(tnames, bar)
  }
}
names(tdata) <- tnames
pdata <- list()
for (d in names(tdata)){
  print(d)
  x <- list(as.data.frame(tdata[[d]]))
  pdata <- append(pdata, x)
  #pdata[[d]] <- as.data.frame(tdata$d)
}
names(pdata) <- tnames
tjson <- toJSON(pdata)
write(tjson, "/gstore/project/hr_brca_heterogeneity/T47D_trace_Seq_v2/Nathan/tran_9545_pcares0.25_BREAST75.json")



##########################
# new cut off for 
# PATHWAYS
# find the biggest gap between normalized absolute change >= ?
# keep all cells in the biggest gap, remove all cells outside boundary

cellbounds <- function(change, scebar){
  norm <- minMax(change)
  position <- which(norm >= 0.15)
  half <- length(change)/2
  poslow <- max(position[position <= half])
  poshigh <- min(position[position > half])
  scebar <- scebar[poslow:poshigh, ]
  return(scebar)
}

bars <- c(bar, bar2, bar3, bar4, bar5, bar6, bar7, bar8, bar9, bar10, bar11)
for (b in bars){
  scetime_bar <- scetime[scetime$barcode == bar11, ]
  changerate <- c()
  for (i in 1:length(scetime_bar$barsurvival)-1){
    changerate <- append(changerate, abs(scetime_bar$barsurvival[i+1] - scetime_bar$barsurvival[i]))
  }
  scebar <- cellbounds(changerate, scetime_bar)
  
  
}


treats <- c('GDC-9545 200nm', 'Baseline')
tdata <- list()
tnames <- c()
##### turn the order by to ptime only or realtime + ptime
for (bar in barlist){
  bar_sce <- sce[, colData(sce)$gfp_bc == bar]
  if (dim(bar_sce)[2] > 25){   #### this line is the minimum number of cells 
    print(bar)
    print(dim(bar_sce)[2])
    df <- colData(bar_sce)[24:length(colData(bar_sce))]
    df <- df[ ,grepl("Score",colnames(df))]
    x <- cbind(bar_sce$slingPseudotime_Avg2, cbind(bar_sce$time_point, bar_sce$barsurvival))
    df <- cbind(x, df)
    colnames(df)[1:3] <- c('Pseudotime', 'RealTime', 'Survival')
    #df <- df[order( df[,2], df[,1] ),] # order by realtime and then order internal by ptime
    df <- df[order(df[,1]), ] # order by pseudotime
    print(dim(df))
    b <- list(df)
    tdata <- append(tdata, b)
    tnames <- append(tnames, bar)
  }
}
names(tdata) <- tnames
pdata <- list()
for (d in names(tdata)){
  print(d)
  x <- list(as.data.frame(tdata[[d]]))
  pdata <- append(pdata, x)
  #pdata[[d]] <- as.data.frame(tdata$d)
}
names(pdata) <- tnames
tjson <- toJSON(pdata)
write(tjson, "/gstore/project/hr_brca_heterogeneity/T47D_trace_Seq_v2/Nathan/tran_9545_pcares0.25_BREAST75.json")
