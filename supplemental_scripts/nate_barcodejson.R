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

# gene names
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
# filter so that at least 25% of cells have 3 counts of included genes
sce <- sce[rowSums(counts(sce) > 1) > ncol(sce)/10, ]
#sce <- sce[rowSums(counts(sce) > 1) > ncol(sce)/10, ]
dim(sce)

###### pca trajectory ######
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
sce$slingPseudotime_Avg2 <- rowMeans(slingPseudotime(sce), na.rm = TRUE)

### ridge plot
# as.character vs not is different view
scetime <- as.character(sce$time_point)  #scetime <- as.character(sce$time_point)
scetime <- cbind(scetime, sce$slingPseudotime_Avg2)  #scetime <- cbind(scetime, as.numeric(sce$slingPseudotime_Avg2))
scetime <- cbind(scetime, as.character(sce$gfp_bc))
colnames(scetime) <- c('time_point', 'pseudotime', 'barcode')
scetime <- as.data.frame(scetime)
scetime$pseudotime <- as.numeric(scetime$pseudotime)

timelab <- c('Day 0', 'Day 1', 'Day 4', 'Day 8', 'Day 26')
ggplot(scetime, aes(x = pseudotime, y = time_point, 
                    group = time_point)) + geom_density_ridges() + 
  theme_ridges() + theme(legend.position = "none") 

dittoPlot(sce, "slingPseudotime_Avg2", group.by = "time_point",
          plots = c("ridgeplot"), ridgeplot.lineweight = 0.5,
          legend.show = FALSE)

# order numerically not by character!!!!
scetime <- scetime[order(scetime[,2]), ]
barsurv <- c()
for (i in 1:length(scetime$barcode)){
  s <- length(which(scetime[c(i:dim(scetime)[1]),]$barcode == scetime$barcode[i]))
  stot <- s / length(scetime[c(i:dim(scetime)[1]),]$barcode)
  print(length(scetime[c(i:dim(scetime)[1]),]$barcode))
  print(s)
  print(stot)
  barsurv <- append(barsurv, stot)
}

scetime$barsurvival <- barsurv
# match order of sce object according to barsurv
scetime2 <- scetime[match(rownames(colData(sce)), rownames(scetime)),]
sce$barsurvival <- scetime2$barsurvival

##################################################
# output json with {Barcode_1: df, Barcode_n: df_n} structure
bar26day <- function(bar, treat_var, time_var){
  x <- barcode_survival[[bar]]
  y <- x[x$treatment == treat_var, ]
  as.numeric(y[y$time_point == time_var, 6])
}


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
    surv <- c()
    for (t in bar_sce$time_point){
      if(t != 'Day 0'){
        surv <- append(surv, bar26day(bar, treats[1], t))
      }
      else {
        surv <- append(surv, bar26day(bar, treats[2], t))
      }
    }
    df <- t(as.data.frame(logcounts(bar_sce)))
    df <- as.data.frame(df)
    x <- cbind(bar_sce$slingPseudotime_Avg2, cbind(bar_sce$time_point, surv))
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
write(tjson, "/gstore/project/hr_brca_heterogeneity/T47D_trace_Seq_v2/Nathan/tran_9545_barcodes4.json")


##################################################
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
    surv <- c()
    for (s in bar_sce$barsurvival){
      surv <- append(surv, s)
    }
    df <- t(as.data.frame(logcounts(bar_sce)))
    df <- as.data.frame(df)
    x <- cbind(bar_sce$slingPseudotime_Avg2, cbind(bar_sce$time_point, surv))
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
write(tjson, "/gstore/project/hr_brca_heterogeneity/T47D_trace_Seq_v2/Nathan/tran_9545_barcodes5.json")

