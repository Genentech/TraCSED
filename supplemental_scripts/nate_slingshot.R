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
library(tradeSeq)
library(biomaRt)
library(org.Hs.eg.db)
library(fgsea)
library(genomitory)
library(RColorBrewer)
library(velociraptor)

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


##############################################
# treatment for trajectory
dim(sce)
treats <- c('Baseline', 'Palbociclib 200nm')
sce <- sce[ , colData(sce)$treatment %in% treats]
dim(sce)
sce <- sce[ , colData(sce)$gfp_bc %in% names(barcode_survival)]
dim(sce)

# do the following for tradeseq at least
# preprocess by gene count
sce <- sce[rowSums(counts(sce) > 3) > ncol(sce)/4, ]
dim(sce)

################################################
# these lines are for specifying a barcode subset condition
survstats <- c()
for (bar in sce$gfp_bc){
  survstats <- append(survstats, bar26day(bar, "Palbociclib 200nm", 'Day 26'))
}
length(survstats)

colData(sce)$bar_surv <- survstats
colData(sce)$bar_label <- ifelse(survstats >= 1, "Growth", "Depleted")
# remove na
sce <- sce[ , !is.na(sce$bar_label)]
dim(sce)

# subset
#sce <- sce[ , colData(sce)$bar_label == 'Growth']
#dim(sce)

#############################################
###### pca ######
set.seed(seed_value)
sce <- runUMAP(sce, dimred = "PCA", BPPARAM = BPPARAM)
dittoDimPlot(sce, "treatment", reduction.use = "UMAP", size = point_size)
dittoDimPlot(sce, "time_point", reduction.use = "UMAP", size = point_size)

# clustering
cluster_resolutions <- c(0.2, 0.3, 0.5, 0.75, 1, 1.25)
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


# cluster data to inform trajectory
clust.data <- as.data.frame(colData(sce)$resolution1)
clust.data <- cbind(clust.data, colData(sce)$time_point)
clust.data <- cbind(clust.data, colData(sce)$bar_label)
colnames(clust.data) <- c('Cluster', 'Time', 'Survival')
table(clust.data)
df <- as.data.frame(table(clust.data))
growth_df <- df[df$Survival == 'Growth', ]
depl_df <- df[df$Survival == 'Depleted', ]
growth_df2 <- reshape(growth_df[,-3], timevar = "Time", idvar = "Cluster", direction = "wide")
row.names(growth_df2) <- growth_df2[,1]
growth_df2 <- growth_df2[,-1]
depl_df2 <- reshape(depl_df[,-3], timevar = "Time", idvar = "Cluster", direction = "wide")
row.names(depl_df2) <- depl_df2[,1]
depl_df2 <- depl_df2[,-1]

ComplexHeatmap::Heatmap(scale(growth_df2), name = 'Scaled No. Cells', 
                        row_names_gp = grid::gpar(fontsize = 4), 
                        column_names_gp = grid::gpar(fontsize = 9),
                        cluster_columns = FALSE)
ComplexHeatmap::Heatmap(scale(depl_df2), name = 'Scaled No. Cells', 
                        row_names_gp = grid::gpar(fontsize = 4), 
                        column_names_gp = grid::gpar(fontsize = 9),
                        cluster_columns = FALSE)


# run slingshot
start_cluster <- "1"
end_cluster <- "12"
sce <- slingshot(sce, reducedDim = "PCA", stretch = 0, approx_points = 200)
sce$slingPseudotime_Global <- sce$slingPseudotime_1
sce$slingPseudotime_1 <- NULL
sce$slingClusters <- NULL


sce <- slingshot(sce, reducedDim = "PCA",
                 clusterLabels = sce$resolution1,
                 start.clus = start_cluster,
                 end.clus = end_cluster,
                 stretch = 0, approx_points = 200)
sce$slingPseudotime_Avg <- rowMeans(slingPseudotime(sce), na.rm = TRUE)

slingLineages(sce)
################################
# UMAPs
# Global
dittoDimPlot(sce, "slingPseudotime_Global", reduction.use = "UMAP",
             size = point_size) + scale_color_viridis_c(direction = -1)
# Avg
dittoDimPlot(sce, "slingPseudotime_Avg", reduction.use = "UMAP",
             size = point_size) + scale_color_viridis_c(direction = -1)
# trajectory
dittoDimPlot(sce, "resolution1", reduction.use = "UMAP", size = point_size,
             do.label = TRUE, legend.show = FALSE,
             add.trajectory.lineages = slingLineages(sce),
             trajectory.cluster.meta = "resolution1")
# growth
dittoDimPlot(sce, "bar_label", reduction.use = "UMAP", size = point_size, 
             legend.show = TRUE)
dittoDimPlot(sce, "bar_surv", reduction.use = "UMAP", size = point_size, 
             legend.show = TRUE)

####################################################
# tradeseq 
sce2 <- sce

lineages <- getLineages(data = reducedDim(sce, type = "UMAP"), clusterLabels = as.numeric(sce$resolution1), 
                        start.clus = start_cluster, end.clus = end_cluster) #define where to start the trajectories
curves <- getCurves(lineages, approx_points = 300, thresh = 0.01, stretch = 0.8, allow.breaks = FALSE, shrink = 0.99)

branch_assignments <- slingBranchID(sce)
pseudotime <- slingPseudotime(curves, na = FALSE)
cellWeights <- slingCurveWeights(curves)
# fit negative binomial GAM
sce <- tradeSeq::fitGAM(as.matrix(counts(sce)), sds = curves, pseudotime = pseudotime, 
                        cellWeights = cellWeights, BPPARAM = BPPARAM)


colData(sce2)$L5_weights <- cellWeights$Lineage5
dittoDimPlot(sce2, "L5_weights", reduction.use = "UMAP", size = point_size, 
             legend.show = TRUE)

# gene names
ids <- row.names(sce)
annots <- select(org.Hs.eg.db, keys=ids, 
                 columns="SYMBOL", keytype="ENSEMBL")
annots <- annots[!duplicated(annots$ENSEMBL), ]
row.names(annots) <- annots$ENSEMBL
annots$GENE <- ifelse(is.na(annots$SYMBOL), annots$ENSEMBL, annots$SYMBOL)
rownames(sce) <- annots$GENE

##########################################
# association test for dynamic expression
ATres <- associationTest(sce)
View(ATres)

topgenes <- rownames(ATres[order(ATres$pvalue), ])[1:20]
pst.ord <- order(sce2$slingPseudotime_1, na.last = NA)
heatdata <- assays(sce)$counts[topgenes, pst.ord]
heatclus <- sce2$resolution1[pst.ord]

heatmap(log1p(heatdata), Colv = NA,
        ColSideColors = brewer.pal(9,"Set1")[heatclus])

unique(brewer.pal(9,"Set1")[heatclus])
unique(heatclus)

plotSmoothers(sce, as.matrix(counts(sce)), "CRIP1")
plotGeneCount(curves, as.matrix(counts(sce)), gene = "CRIP1")
plotGeneCount(curves, as.matrix(counts(sce)), clusters = sce$resolution0.75, models = sce)

##########################################
# start vs end of lineage
startRes <- startVsEndTest(sce)
View(startRes)

topgenes <- rownames(startRes[order(startRes$waldStat, decreasing=TRUE), ])[1:20]
pst.ord <- order(sce2$slingPseudotime_1, na.last = NA)
heatdata <- assays(sce)$counts[topgenes, pst.ord]
heatclus <- sce2$resolution0.75[pst.ord]
heatmap(log1p(heatdata), Colv = NA,
        ColSideColors = brewer.pal(9,"Set1")[heatclus])

unique(brewer.pal(9,"Set1")[heatclus])
unique(heatclus)

plotSmoothers(sce, as.matrix(counts(sce)), "COX1")
plotGeneCount(curves, as.matrix(counts(sce)), gene = "COX1")

##########################################
# are end states different between lineages
endRes <- diffEndTest(sce)
View(endRes)

topgenes <- rownames(endRes[order(endRes$waldStat, decreasing=TRUE), ])[1:20]
pst.ord <- order(sce2$slingPseudotime_1, na.last = NA)
heatdata <- assays(sce)$counts[topgenes, pst.ord]
heatclus <- sce2$resolution0.75[pst.ord]
heatmap(log1p(heatdata), Colv = NA,
        ColSideColors = brewer.pal(9,"Set1")[heatclus])

unique(brewer.pal(9,"Set1")[heatclus])
unique(heatclus)

plotSmoothers(sce, as.matrix(counts(sce)), "PTMA")
plotGeneCount(curves, as.matrix(counts(sce)), gene = "PTMA")



dittoDimPlot(sce2, "slingPseudotime_1", reduction.use = "UMAP",
             size = point_size) + scale_color_viridis_c(direction = -1)



####################################################
# Run scVelo with steady state
set.seed(seed_value)
sce_velo_steady_state <- scvelo(sce,
                                mode = "steady_state",
                                assay.X = "counts",
                                assay.spliced = "counts_spliced",
                                assay.unspliced = "counts_unspliced",
                                use.dimred = "PCA",
                                BPPARAM = BPPARAM)
# Add RNA velocity information to the SingleCellExperiment object
sce$velocity_pseudotime_steady_state <- sce_velo_steady_state$velocity_pseudotime
sce$velocity_length_steady_state <- sce_velo_steady_state$velocity_length
sce$velocity_confidence_steady_state <- sce_velo_steady_state$velocity_confidence
# Project velocity vectors to the UMAP coordinates
embedded_steady_state <- embedVelocity(reducedDim(sce, "UMAP"),
                                       sce_velo_steady_state)
grid_df_steady_state <- gridVectors(reducedDim(sce, "UMAP"),
                                    embedded_steady_state, resolution = 30)


dittoDimPlot(sce, "velocity_pseudotime_steady_state",
             reduction.use = "UMAP", size = point_size) +
  geom_segment(
    data = grid_df_steady_state,
    mapping = aes(x = start.1, y = start.2, xend = end.1, yend = end.2),
    arrow = arrow(length = unit(0.05, "inches"), type = "open")
  )

dittoDimPlot(sce, "velocity_pseudotime_steady_state",
             reduction.use = "UMAP", size = point_size)

dittoPlot(sce, "velocity_pseudotime_steady_state", group.by = "resolution0.75",
          plots = c("ridgeplot"), ridgeplot.lineweight = 0.5,
          legend.show = FALSE)
dittoDimPlot(sce, "velocity_confidence_steady_state",
             reduction.use = "UMAP", size = point_size)



# Run scVelo with dynamical mode
set.seed(seed_value)
sce_velo_dynamical <- scvelo(sce,
                             mode = "dynamical",
                             assay.X = "counts",
                             assay.spliced = "counts_spliced",
                             assay.unspliced = "counts_unspliced",
                             use.dimred = "PCA",
                             BPPARAM = BPPARAM)
# Add RNA velocity information to the SingleCellExperiment object
sce$velocity_pseudotime_dynamical <- sce_velo_dynamical$velocity_pseudotime
sce$velocity_length_dynamical <- sce_velo_dynamical$velocity_length
sce$velocity_confidence_dynamical <- sce_velo_dynamical$velocity_confidence
# Project velocity vectors to the UMAP coordinates
embedded_dynamical <- embedVelocity(reducedDim(sce, "UMAP"),
                                    sce_velo_dynamical)
grid_df_dynamical <- gridVectors(reducedDim(sce, "UMAP"),
                                 embedded_dynamical, resolution = 30)

dittoDimPlot(sce, "velocity_pseudotime_dynamical",
             reduction.use = "UMAP", size = point_size) +
  geom_segment(
    data = grid_df_dynamical,
    mapping = aes(x = start.1, y = start.2, xend = end.1, yend = end.2),
    arrow = arrow(length = unit(0.05, "inches"), type = "open")
  )

dittoDimPlot(sce, "velocity_pseudotime_dynamical",
             reduction.use = "UMAP", size = point_size)

dittoPlot(sce, "velocity_pseudotime_dynamical", group.by = "resolution0.75",
          plots = c("ridgeplot"), ridgeplot.lineweight = 0.5,
          legend.show = FALSE)
dittoDimPlot(sce, "velocity_confidence_dynamical",
             reduction.use = "UMAP", size = point_size)

# gene names
ids <- row.names(sce)
annots <- select(org.Hs.eg.db, keys=ids, 
                 columns="SYMBOL", keytype="ENSEMBL")
annots <- annots[!duplicated(annots$ENSEMBL), ]
row.names(annots) <- annots$ENSEMBL
annots$GENE <- ifelse(is.na(annots$SYMBOL), annots$ENSEMBL, annots$SYMBOL)
rownames(sce) <- annots$GENE

x <- gene("KIF1B", object = sce, assay = "counts_spliced")
y <- gene("KIF1B", object = sce, assay = "counts_unspliced")
plot(x, y, pch = 19, col = "black", xlab = 'Spliced', ylab = 'Unspliced')

###################################################
###### plsr prep ######
###################################################
# filter for plsr model
# PLSR
set.seed(1)
######### filter the data
dim(sce)
# only cells in our survival list
sce <- sce[ ,colData(sce)$gfp_bc %in% names(barcode_survival)]
dim(sce)
# filter so that at least 25% of cells have 3 counts of included genes
sce <- sce[rowSums(counts(sce) > 3) > ncol(sce)/4, ]
dim(sce)


# model data
df2 <- t(as.data.frame(logcounts(sce)))
df2 <- as.data.frame(df2)
r <- cbind(day1surv2, cbind(day4surv2, cbind(day8surv2, day26surv2)))
x <- as.vector(colData(sce)$time_point)
x <- replace(x, x == "Day 0", 0)
x <- replace(x, x == 'Day 1', 1)
x <- replace(x, x == 'Day 4', 4)
x <- replace(x, x == 'Day 8', 8)
x <- replace(x, x == 'Day 26', 26)
#x <- as.numeric(x)

model <- plsr(x ~ ., data=df2, ncomp=30, scale=TRUE, validation='CV')
model <- plsr(r ~ ., data=df2, ncomp=30, scale=TRUE, validation='CV')
summary(model)
plot(model, plottype = "scores", comps = 1:5)


# add plsr components
dim(reducedDim(sce, "PCA"))
plsr_scores <- scores(model)[,1:5] ### change slice to get more components
reducedDim(sce, "PLSR") <- plsr_scores
dim(reducedDim(sce, "PLSR"))

# run umap
set.seed(1)
sce <- runUMAP(sce, dimred = "PLSR", BPPARAM = BPPARAM)

# treatment umap
dittoDimPlot(sce, "treatment", reduction.use = "UMAP", size = point_size)
dittoDimPlot(sce, "time_point", reduction.use = "UMAP", size = point_size)

# clustering
cluster_resolutions <- c(0.15, 0.2, 0.3, 0.5, 0.75, 1, 1.25)
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

##############################
# run slingshot
start_cluster <- "2"
end_cluster <- "6"
sce <- slingshot(sce, reducedDim = "PLSR", stretch = 0, approx_points = 200)
sce$slingPseudotime_Global <- sce$slingPseudotime_1
sce$slingPseudotime_1 <- NULL
sce$slingClusters <- NULL


sce <- slingshot(sce, reducedDim = "PLSR",
                 clusterLabels = sce$resolution0.3,
                 start.clus = start_cluster,
                 end.clus = end_cluster,
                 stretch = 0, approx_points = 200)
sce$slingPseudotime_Avg <- rowMeans(slingPseudotime(sce), na.rm = TRUE)

slingLineages(sce)
################################
# UMAPs
# Global
dittoDimPlot(sce, "slingPseudotime_Global", reduction.use = "UMAP",
             size = point_size) + scale_color_viridis_c(direction = -1)
# Avg
dittoDimPlot(sce, "slingPseudotime_Avg", reduction.use = "UMAP",
             size = point_size) + scale_color_viridis_c(direction = -1)
# trajectory
dittoDimPlot(sce, "resolution0.3", reduction.use = "UMAP", size = point_size,
             do.label = TRUE, legend.show = FALSE,
             add.trajectory.lineages = slingLineages(sce),
             trajectory.cluster.meta = "resolution0.3")

lineages <- getLineages(data = reducedDim(sce, type = "UMAP"), clusterLabels = as.numeric(sce$resolution0.3), 
                        start.clus = start_cluster, end.clus = end_cluster) #define where to start the trajectories
curves <- getCurves(lineages, approx_points = 300, thresh = 0.01, stretch = 0.8, allow.breaks = FALSE, shrink = 0.99)

branch_assignments <- slingBranchID(sce)
pseudotime <- slingPseudotime(curves, na = FALSE)
cellWeights <- slingCurveWeights(curves)
cellWeights <- as.data.frame(cellWeights)

colData(sce)$L2_weights <- cellWeights$Lineage2
dittoDimPlot(sce, "L2_weights", reduction.use = "UMAP", size = point_size, 
             legend.show = TRUE)
plotGeneCount(curves, as.matrix(counts(sce)), clusters = sce$resolution0.3, models = sce)

# growth
dittoDimPlot(sce, "bar_label", reduction.use = "UMAP", size = point_size, 
             legend.show = TRUE)
dittoDimPlot(sce, "bar_surv", reduction.use = "UMAP", size = point_size, 
             legend.show = TRUE)


# fit negative binomial GAM
sce <- tradeSeq::fitGAM(as.matrix(counts(sce)), sds = curves, pseudotime = pseudotime, 
                        cellWeights = cellWeights, BPPARAM = BPPARAM)


# gene names
ids <- row.names(sce)
annots <- select(org.Hs.eg.db, keys=ids, 
                 columns="SYMBOL", keytype="ENSEMBL")
annots <- annots[!duplicated(annots$ENSEMBL), ]
row.names(annots) <- annots$ENSEMBL
annots$GENE <- ifelse(is.na(annots$SYMBOL), annots$ENSEMBL, annots$SYMBOL)
rownames(sce) <- annots$GENE

##########################################
# association test for dynamic expression
ATres <- associationTest(sce)
View(ATres)

topgenes <- rownames(ATres[order(ATres$pvalue), ])[1:25]
pst.ord <- order(sce$slingPseudotime_1, na.last = NA)
heatdata <- assays(sce)$counts[topgenes, pst.ord]
heatclus <- sce$resolution0.3[pst.ord]

heatmap(log1p(heatdata), Colv = NA,
        ColSideColors = brewer.pal(9,"Set1")[heatclus])

unique(brewer.pal(9,"Set1")[heatclus])
unique(heatclus)

plotSmoothers(sce, as.matrix(counts(sce)), "NCAPD3")
plotGeneCount(curves, as.matrix(counts(sce)), gene = "NCAPD3")
plotGeneCount(curves, as.matrix(counts(sce)), clusters = sce$resolution0.75, models = sce)

##########################################
# start vs end of lineage
startRes <- startVsEndTest(sce)
View(startRes)

topgenes <- rownames(startRes[order(startRes$waldStat, decreasing=TRUE), ])[1:20]
pst.ord <- order(sce$slingPseudotime_1, na.last = NA)
heatdata <- assays(sce)$counts[topgenes, pst.ord]
heatclus <- sce$resolution0.3[pst.ord]
heatmap(log1p(heatdata), Colv = NA,
        ColSideColors = brewer.pal(9,"Set1")[heatclus])

unique(brewer.pal(9,"Set1")[heatclus])
unique(heatclus)

plotSmoothers(sce, as.matrix(counts(sce)), "IL11RA")
plotGeneCount(curves, as.matrix(counts(sce)), gene = "IL11RA")

##########################################
# are end states different between lineages
endRes <- diffEndTest(sce)
View(endRes)

topgenes <- rownames(endRes[order(endRes$waldStat, decreasing=TRUE), ])[1:20]
pst.ord <- order(sce2$slingPseudotime_1, na.last = NA)
heatdata <- assays(sce)$counts[topgenes, pst.ord]
heatclus <- sce2$resolution0.75[pst.ord]
heatmap(log1p(heatdata), Colv = NA,
        ColSideColors = brewer.pal(9,"Set1")[heatclus])

unique(brewer.pal(9,"Set1")[heatclus])
unique(heatclus)

plotSmoothers(sce, as.matrix(counts(sce)), "H2AZ1")
plotGeneCount(curves, as.matrix(counts(sce)), gene = "H2AZ1")



dittoDimPlot(sce2, "slingPseudotime_1", reduction.use = "UMAP",
             size = point_size) + scale_color_viridis_c(direction = -1)
