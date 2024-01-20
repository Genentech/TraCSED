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
library(DESeq2)
library(biomaRt)
library(org.Hs.eg.db)
library(fgsea)
library(genomitory)

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

# pathways of interest
availcol <- as.data.frame(availableFiles(type="collection"))
#View(availcol)
pid_pathways <- getFeatureSetCollection("GMTY42:human/C2.CP.PID.gmt.bz2@REVISION-1")
cgp_pathways <- getFeatureSetCollection("GMTY42:human/C2.CGP.gmt.bz2@REVISION-1")
c2_pathways <- getFeatureSetCollection("GMTY42:human/C2.gmt.bz2@REVISION-1")
breast_er_pathways <- getFeatureSetCollection("GMTY194:analysis/breast.gmt.bz2@REVISION-2")

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

##############################################################################
# attempt 1 - day 0 surviving resistant vs surviving depleted
assayNames(sce)
dim(logcounts(sce))
sce2 <- sce[ , colData(sce)$time_point == 'Day 0']
dim(logcounts(sce2))
sce2 <- sce2[ ,colData(sce2)$gfp_bc %in% barcode_summary$gfp_bc]
dim(logcounts(sce2))

## Remove lowly expressed genes which have less than 10 cells with any counts
sce2 <- sce2[rowSums(counts(sce2) > 1) >= 10, ]
dim(logcounts(sce2))
exp <- aggregate(t(counts(sce2)), by=list(colData(sce2)[,c('gfp_bc')]), FUN=sum)

###### get survival conditions

bars <- as.data.frame(exp$Group.1)

bar26day <- function(bar, treat_var){
  x <- barcode_survival[[bar]]
  y <- x[x$treatment == treat_var, ]
  as.numeric(y[y$time_point == 'Day 26', 6])
}
survstats <- c()
for (bar in exp$Group.1){
  if (bar %in% names(barcode_survival)){
    print(bar26day(bar, "GDC-9545 + Palbociclib"))
    survstats <- append(survstats, bar26day(bar, "GDC-9545 + Palbociclib"))
  }
  else {
    print(bar)
    print('killed')
    survstats <- append(survstats, -Inf)
  }
}

condition <- cbind(bars, ifelse(survstats >= 0, "Growth", ifelse(survstats == -Inf, 'Remove', 'Depletion')))
row.names(exp) <- bars[,1]
exp <- exp [, -1]
exp <- t(exp)
row.names(condition) <- bars[,1]
colnames(condition) <- c('barcode', 'condition')

# prune out barcodes with "remove" label
condition <- condition[!condition$condition == 'Remove', ]
exp <- exp[,colnames(exp) %in% row.names(condition)]



# Deseq data
dds <- DESeqDataSetFromMatrix(exp, condition, ~ condition)
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "Depletion", "Growth"))
res <- res[order(res$padj), ]
resdata <- as.data.frame(res)
ids <- row.names(resdata)
annots <- select(org.Hs.eg.db, keys=ids, 
                 columns="SYMBOL", keytype="ENSEMBL")
annots <- annots[!duplicated(annots$ENSEMBL), ]
row.names(annots) <- annots$ENSEMBL
resdata <- cbind(annots$SYMBOL, resdata)
colnames(resdata)[1] <- 'Gene Symbol'
siggenes <- resdata[resdata$padj <= 0.05, ]
siggenes <- siggenes[complete.cases(siggenes), ]
plotCounts(dds, 'ENSG00000165813', intgroup = "condition")




##############################################################################
# attempt 2 - day 0 surviving resistant vs surviving and killed depleted
assayNames(sce)
dim(logcounts(sce))
sce <- sce[ ,colData(sce)$gfp_bc %in% names(barcode_survival)]
dim(logcounts(sce))
## Remove lowly expressed genes which have less than 10 cells with any counts
sce <- sce[rowSums(counts(sce) > 1) >= 10, ]
dim(logcounts(sce))
sce <- sce[ ,colData(sce)$time_point == 'Day 0']
dim(sce)
groups <- names(barcode_survival)
exp <- exp <- aggregate(t(counts(sce)), by=list(colData(sce)[,c('gfp_bc')]), FUN=sum)
###### get survival conditions
# in this case binarize the data by if by day 26 
# the combination had more or less cells that day 0
bar26day <- function(bar, treat_var){
  x <- barcode_survival[[bar]]
  y <- x[x$treatment == treat_var, ]
  as.numeric(y[y$time_point == 'Day 26', 6])
}
survstats <- c()
for (bar in exp$Group.1){
  survstats <- append(survstats, bar26day(bar, "GDC-9545 200nm"))
}
bars <- as.data.frame(exp$Group.1)
condition <- cbind(bars, ifelse(survstats >= 0, "Growth", "Depletion"))
row.names(exp) <- bars[,1]
exp <- exp [, -1]
exp <- t(exp)

condition <- as.data.frame(ifelse(survstats >= 0, "Growth", "Depletion"))
row.names(condition) <- bars[,1]
colnames(condition) <- 'condition'
# Deseq data
dds <- DESeqDataSetFromMatrix(exp, condition, ~ condition)
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "Depletion", "Growth"))
res <- res[order(res$padj), ]
resdata <- as.data.frame(res)
ids <- row.names(resdata)
annots <- select(org.Hs.eg.db, keys=ids, 
                 columns="SYMBOL", keytype="ENSEMBL")
annots <- annots[!duplicated(annots$ENSEMBL), ]
row.names(annots) <- annots$ENSEMBL
resdata <- cbind(annots$SYMBOL, resdata)
colnames(resdata)[1] <- 'Gene Symbol'
plotCounts(dds, 'ENSG00000159763', intgroup = "condition", main = 'PIP - GDC-9545')
# ENSG00000163347 - CLDN1
# ENSG00000266402 - SNHG25
# ENSG00000173267 - SNCG

##############################################################################
# attempt 3 - day 0 surviving resistant vs killed depleted
assayNames(sce)
dim(logcounts(sce))
sce2 <- sce[ , colData(sce)$time_point == 'Day 0']
dim(logcounts(sce2))
sce2 <- sce2[ ,colData(sce2)$gfp_bc %in% barcode_summary$gfp_bc]
dim(logcounts(sce2))

## Remove lowly expressed genes which have less than 10 cells with any counts
sce2 <- sce2[rowSums(counts(sce2) > 1) >= 10, ]
dim(logcounts(sce2))
exp <- aggregate(t(counts(sce2)), by=list(colData(sce2)[,c('gfp_bc')]), FUN=sum)

###### get survival conditions
# in this case binarize the data by if by day 26 
# the combination had more or less cells that day 0

# if in barcode_summary day26 entropy > 0
bars <- as.data.frame(exp$Group.1)

bar26day <- function(bar, treat_var){
  x <- barcode_survival[[bar]]
  y <- x[x$treatment == treat_var, ]
  as.numeric(y[y$time_point == 'Day 26', 6])
}
survstats <- c()
for (bar in exp$Group.1){
  if (bar %in% names(barcode_survival)){
    print(bar26day(bar, "GDC-9545 + Palbociclib"))
    survstats <- append(survstats, bar26day(bar, "GDC-9545 + Palbociclib"))
  }
  else {
    print(bar)
    print('killed')
    survstats <- append(survstats, -Inf)
  }
}

condition <- cbind(bars, ifelse(survstats >= 0, "Growth", ifelse(survstats == -Inf, 'Depletion', 'Remove')))
row.names(exp) <- bars[,1]
exp <- exp [, -1]
exp <- t(exp)
row.names(condition) <- bars[,1]
colnames(condition) <- c('barcode', 'condition')

# prune out barcodes with "remove" label
condition <- condition[!condition$condition == 'Remove', ]
exp <- exp[,colnames(exp) %in% row.names(condition)]

# Deseq data
dds <- DESeqDataSetFromMatrix(exp, condition, ~ condition)
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "Depletion", "Growth"))
res <- res[order(res$padj), ]
resdata <- as.data.frame(res)
ids <- row.names(resdata)
annots <- select(org.Hs.eg.db, keys=ids, 
                 columns="SYMBOL", keytype="ENSEMBL")
annots <- annots[!duplicated(annots$ENSEMBL), ]
row.names(annots) <- annots$ENSEMBL
resdata <- cbind(annots$SYMBOL, resdata)
colnames(resdata)[1] <- 'Gene Symbol'
siggenes <- resdata[resdata$padj <= 0.05, ]
siggenes <- siggenes[complete.cases(siggenes), ]
plotCounts(dds, 'ENSG00000165813', intgroup = "condition")









##############################################################################
# attempt 4 - attempt 3, but added cells at day0 filter and gsea
trx <- "Palbociclib 200nm"
itimevar <- 'Day 0'
endvar <- 'Day 26'
assayNames(sce)
dim(logcounts(sce))

day0count <- c()
for (bar in barcode_summary$gfp_bc){
  x <- barcode_survival[[bar]]
  day0count <- append(day0count, as.numeric(x[x$time_point == itimevar, 3]))
}
barcode_summary$Day0cellcount <- day0count
barcode_summary2 <- barcode_summary[barcode_summary$Day0cellcount >= 5, ]

sce2 <- sce[ , colData(sce)$time_point == itimevar]
dim(logcounts(sce2))
sce2 <- sce2[ ,colData(sce2)$gfp_bc %in% barcode_summary2$gfp_bc]
dim(logcounts(sce2))

## Remove lowly expressed genes which have less than 10 cells with any counts
sce2 <- sce2[rowSums(counts(sce2) > 1) >= 10, ]
dim(logcounts(sce2))
exp <- aggregate(t(counts(sce2)), by=list(colData(sce2)[,c('gfp_bc')]), FUN=sum)

###### get survival conditions
# in this case binarize the data by if by day 26 
# the combination had more or less cells that day 0

bars <- as.data.frame(exp$Group.1)

bar26day <- function(bar, treat_var){
  x <- barcode_survival[[bar]]
  y <- x[x$treatment == treat_var, ]
  as.numeric(y[y$time_point == endvar, 6])
}
survstats <- c()
for (bar in exp$Group.1){
  if (bar %in% names(barcode_survival)){
    survstats <- append(survstats, bar26day(bar, trx))
  }
  else {
    print(bar)
    print('killed')
    survstats <- append(survstats, -Inf)
  }
}

condition <- cbind(bars, ifelse(survstats >= 0, "Growth", ifelse(survstats == -Inf, 'Depletion', 'Remove')))
row.names(exp) <- bars[,1]
exp <- exp [, -1]
exp <- t(exp)
row.names(condition) <- bars[,1]
colnames(condition) <- c('barcode', 'condition')

# prune out barcodes with "remove" label
condition <- condition[!condition$condition == 'Remove', ]
exp <- exp[,colnames(exp) %in% row.names(condition)]

# Deseq data
dds <- DESeqDataSetFromMatrix(exp, condition, ~ condition)
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "Depletion", "Growth"))
res <- res[order(res$padj), ]
resdata <- as.data.frame(res)
ids <- row.names(resdata)
annots <- select(org.Hs.eg.db, keys=ids, 
                 columns="SYMBOL", keytype="ENSEMBL")
annots <- annots[!duplicated(annots$ENSEMBL), ]
row.names(annots) <- annots$ENSEMBL
resdata <- cbind(annots$SYMBOL, resdata)
colnames(resdata)[1] <- 'Gene Symbol'
siggenes <- resdata[resdata$padj <= 0.2, ]
siggenes <- siggenes[complete.cases(siggenes), ]
plotCounts(dds, 'ENSG00000143222', intgroup = "condition", main = 'UFC1')
plotCounts(dds, 'ENSG00000005020', intgroup = "condition", main = 'SKAP2')
plotCounts(dds, 'ENSG00000109321', intgroup = "condition", main = 'AREG')
plotCounts(dds, 'ENSG00000159763', intgroup = "condition", main = 'PIP')


#example run - use of FC vs pval
# see https://bioinformatics-core-shared-training.github.io/cruk-summer-school-2018/RNASeq2018/html/06_Gene_set_testing.nb.html#fgsea

# c2 pathways - get fold change for enrichment
#ranks <- siggenes$log2FoldChange
#ranks <- ifelse(resdata$padj <= 0.2, resdata$log2FoldChange, 0)
#names(ranks) <- row.names(resdata)
ranks <- siggenes$log2FoldChange
names(ranks) <- row.names(siggenes)
head(ranks)
barplot(sort(ranks, decreasing = T))
pathways_genes <- as.list(c2_pathways)
paths <- as.data.frame(elementMetadata(c2_pathways)[1])
names(pathways_genes) <- paths$name
fgsea_c2 <- fgsea(pathways_genes, ranks, minSize = 10, maxSize = 500)

plotEnrichment(pathways_genes[["GOZGIT_ESR1_TARGETS_DN"]], ranks)


# pid pathways - get fold change for enrichment
pathways_genes <- as.list(pid_pathways)
paths <- as.data.frame(elementMetadata(pid_pathways)[1])
names(pathways_genes) <- paths$name
fgsea_pid <- fgsea(pathways_genes, ranks, minSize = 10, maxSize = 500)

plotEnrichment(pathways_genes[["GOZGIT_ESR1_TARGETS_DN"]], ranks)


# cgp pathways - get fold change for enrichment
pathways_genes <- as.list(cgp_pathways)
paths <- as.data.frame(elementMetadata(cgp_pathways)[1])
names(pathways_genes) <- paths$name
fgsea_cgp <- fgsea(pathways_genes, ranks, minSize = 10, maxSize = 500)

plotEnrichment(pathways_genes[["GOZGIT_ESR1_TARGETS_DN"]], ranks)

##############################################################################
# attempt 5 - attempt 2 (all) with cell filter and gsea
trx <- 'GDC-9545 200nm'
itimevar <- 'Day 0'
endvar <- 'Day 26'
assayNames(sce)
dim(logcounts(sce))

day0count <- c()
for (bar in barcode_summary$gfp_bc){
  x <- barcode_survival[[bar]]
  day0count <- append(day0count, as.numeric(x[x$time_point == itimevar, 3]))
}
barcode_summary$Day0cellcount <- day0count
barcode_summary2 <- barcode_summary[barcode_summary$Day0cellcount >= 5, ]

sce2 <- sce[ ,colData(sce)$time_point == itimevar]
dim(sce2)
sce2 <- sce2[ ,colData(sce2)$gfp_bc %in% barcode_summary2$gfp_bc]
dim(logcounts(sce2))
## Remove lowly expressed genes which have less than 10 cells with any counts
sce2 <- sce2[rowSums(counts(sce2) > 1) >= 10, ]
dim(logcounts(sce2))
exp <- aggregate(t(counts(sce2)), by=list(colData(sce2)[,c('gfp_bc')]), FUN=sum)
###### get survival conditions
# in this case binarize the data by if by day 26 
# the combination had more or less cells that day 0
bar26day <- function(bar, treat_var){
  x <- barcode_survival[[bar]]
  y <- x[x$treatment == treat_var, ]
  as.numeric(y[y$time_point == endvar, 6])
}
survstats <- c()
for (bar in exp$Group.1){
  survstats <- append(survstats, bar26day(bar, trx))
}
bars <- as.data.frame(exp$Group.1)
condition <- cbind(bars, ifelse(survstats >= 0, "Growth", "Depletion"))
row.names(exp) <- bars[,1]
exp <- exp [, -1]
exp <- t(exp)

condition <- as.data.frame(ifelse(survstats >= 0, "Growth", "Depletion"))
row.names(condition) <- bars[,1]
colnames(condition) <- 'condition'
# Deseq data
dds <- DESeqDataSetFromMatrix(exp, condition, ~ condition)
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "Depletion", "Growth"))
res <- res[order(res$padj), ]
resdata <- as.data.frame(res)
ids <- row.names(resdata)
annots <- select(org.Hs.eg.db, keys=ids, 
                 columns="SYMBOL", keytype="ENSEMBL")
annots <- annots[!duplicated(annots$ENSEMBL), ]
row.names(annots) <- annots$ENSEMBL
resdata <- cbind(annots$SYMBOL, resdata)
colnames(resdata)[1] <- 'Gene Symbol'
siggenes <- resdata[resdata$padj <= 0.2, ]
siggenes <- siggenes[complete.cases(siggenes), ]
plotCounts(dds, 'ENSG00000266402', intgroup = "condition", main='SNHG25')

# fgsea analysis
ranks <- gseaDat$log2FoldChange
names(ranks) <- gseaDat$`Gene Symbol`
head(ranks)
barplot(sort(ranks, decreasing = T))

# pathway params
pathway_file <- "~/git/HTSEQ-1022-analysis-pipeline-for-trace-seq/data/pathway_genes_NGS4327.tab"
gene_id_column <- "symbol"
point_size <- 0.1
# prep geneset data
pathway_df <- read.delim(pathway_file, header = FALSE)
pathway_genes <- strsplit(pathway_df[, 2], ",")
names(pathway_genes) <- pathway_df[, 1]

##### use FC and make 0 fold change for genes not significant
#example run - use of FC vs pval
# see https://bioinformatics-core-shared-training.github.io/cruk-summer-school-2018/RNASeq2018/html/06_Gene_set_testing.nb.html#fgsea
samp <- pathway_genes$MAPK
samp <- append(samp, pathway_genes$HALLMARK_E2F_TARGETS)
fcranks <- resdata$log2FoldChange[1:length(samp)]
padjranks <- resdata$padj[1:length(samp)]
names(fcranks) <- samp
names(padjranks) <- samp
# run with pathways
fgsea_fc <- fgsea(pathway_genes, fcranks, minSize=10, maxSize = 500)
fgsea_padj <- fgsea(pathway_genes, padjranks, minSize=10, maxSize = 500)

plotEnrichment(pathway_genes[["HALLMARK_E2F_TARGETS"]], fcranks)
plotEnrichment(pathway_genes[["HALLMARK_E2F_TARGETS"]], padjranks)



####################################
# gsea pathways
#example run - use of FC vs pval
# see https://bioinformatics-core-shared-training.github.io/cruk-summer-school-2018/RNASeq2018/html/06_Gene_set_testing.nb.html#fgsea
samp <- pathway_genes$MAPK
samp <- append(samp, pathway_genes$HALLMARK_E2F_TARGETS)
fcranks <- resdata$log2FoldChange[1:length(samp)]
padjranks <- resdata$padj[1:length(samp)]
names(fcranks) <- samp
names(padjranks) <- samp
# run with pathways
fgsea_fc <- fgsea(pathway_genes, fcranks, minSize=10, maxSize = 500)
fgsea_padj <- fgsea(pathway_genes, padjranks, minSize=10, maxSize = 500)

plotEnrichment(pathway_genes[["HALLMARK_E2F_TARGETS"]], fcranks)
plotEnrichment(pathway_genes[["HALLMARK_E2F_TARGETS"]], padjranks)



#########################################################
# pseudotime bulk DE - group at Day X above and below ptime threshold
scei <- sce[ ,colData(sce)$time_point == 'Day 0']
dim(scei)
min(scei$slingPseudotime_Avg2)
max(scei$slingPseudotime_Avg2)
median(scei$slingPseudotime_Avg2)
divline <- ifelse(scei$slingPseudotime_Avg2 >= 40, 'Above', 'Below')
scei$divline <- divline 


markers.div <- scoreMarkers(scei, groups = scei$divline)
above <- as.data.frame(markers.div$Above)
ranks <- above$mean.logFC.cohen
names(ranks) <- row.names(above)

# c2 pathways
c2_pathways <- getFeatureSetCollection("GMTY42:human/C2.gmt.bz2@REVISION-1")
pathways_genes <- as.list(c2_pathways)
paths <- as.data.frame(elementMetadata(c2_pathways)[1])
names(pathways_genes) <- paths$name
fgsea_c2 <- fgsea(pathways_genes, ranks, minSize = 10, maxSize = 500)

# hallmark pathways
hallmark_pathways <- getFeatureSetCollection('GMTY42:human/H.gmt.bz2@REVISION-1')
pathways_genes <- as.list(hallmark_pathways)
paths <- as.data.frame(elementMetadata(hallmark_pathways)[1])
names(pathways_genes) <- paths$name
fgsea_hall <- fgsea(pathways_genes, ranks, minSize = 10, maxSize = 500)


#########################################################
# cellcycle bulk DE - group at Day X in cellcycle stage Y for depleted and growing clones
trx <- 'GDC-9545 200nm'
endvar <- 'Day 26'
phase <- 'G1'
sce <- sce[ ,colData(sce)$treatment == trx]
sce <- sce[ ,colData(sce)$time_point == endvar]
sce <- sce[ ,colData(sce)$phase == phase]
dim(sce)
sce <- sce[rowSums(counts(sce) > 3) >= 10, ]
dim(sce)

###### get survival conditions
# in this case binarize the data by if by day 26 
# the combination had more or less cells that day 0
bar26day <- function(bar, treat_var){
  x <- barcode_survival[[bar]]
  y <- x[x$treatment == treat_var, ]
  as.numeric(y[y$time_point == endvar, 6])
}
survstats <- c()
for (bar in sce$gfp_bc){
  if (bar %in% names(barcode_survival)){
    survstats <- append(survstats, bar26day(bar, trx))
  }
  else {
    survstats <- append(survstats, -Inf)
  }
}
sce$surv <- ifelse(survstats >= 0, 'Growth', 'Depletion')


markers.div <- scoreMarkers(sce, groups = sce$surv)
growth <- as.data.frame(markers.div$Growth)
ranks <- growth$mean.logFC.cohen
names(ranks) <- row.names(growth)

# c2 pathways
c2_pathways <- getFeatureSetCollection("GMTY42:human/C2.gmt.bz2@REVISION-1")
pathways_genes <- as.list(c2_pathways)
paths <- as.data.frame(elementMetadata(c2_pathways)[1])
names(pathways_genes) <- paths$name
fgsea_c2 <- fgsea(pathways_genes, ranks, minSize = 10, maxSize = 500)

# hallmark pathways
hallmark_pathways <- getFeatureSetCollection('GMTY42:human/H.gmt.bz2@REVISION-1')
pathways_genes <- as.list(hallmark_pathways)
paths <- as.data.frame(elementMetadata(hallmark_pathways)[1])
names(pathways_genes) <- paths$name
fgsea_hall <- fgsea(pathways_genes, ranks, minSize = 10, maxSize = 500)

