setwd('/gstore/project/hr_brca_heterogeneity/T47D_trace_Seq_v2/Nathan/barmodels_ptime7')

library(biomaRt)
library(org.Hs.eg.db)
library(fgsea)
library(genomitory)
library(EnsDb.Hsapiens.v79)

pid_pathways <- getFeatureSetCollection("GMTY42:human/C2.CP.PID.gmt.bz2@REVISION-1")
cgp_pathways <- getFeatureSetCollection("GMTY42:human/C2.CGP.gmt.bz2@REVISION-1")
c2_pathways <- getFeatureSetCollection("GMTY42:human/C2.gmt.bz2@REVISION-1")
breast_er_pathways <- getFeatureSetCollection("GMTY194:analysis/breast.gmt.bz2@REVISION-2")

##########################
# transformer
transbar <- read.csv('transformer_GFPBC_libB_37798_feature_importance.csv')

geneIDs <- ensembldb::select(EnsDb.Hsapiens.v79, keys= transbar$feature, keytype = "SYMBOL", columns = c("SYMBOL","GENEID"))
geneIDs <- geneIDs[!duplicated(geneIDs$SYMBOL),]
transbar$geneid <- ifelse(transbar$feature %in% geneIDs$SYMBOL, geneIDs$GENEID, transbar$feature)

transranks <- transbar$fc_mae
names(transranks) <- transbar$geneid

# hallmark pathways
hallmark_pathways <- getFeatureSetCollection('GMTY42:human/H.gmt.bz2@REVISION-1')
pathways_genes <- as.list(hallmark_pathways)
paths <- as.data.frame(elementMetadata(hallmark_pathways)[1])
names(pathways_genes) <- paths$name
trans_fgsea_hall <- fgsea(pathways_genes, transranks, minSize = 10, maxSize = 500)

# percent content of pathway in data
# cgp pathways filtered 
cgp_pathways <- getFeatureSetCollection("GMTY42:human/C2.CGP.gmt.bz2@REVISION-1")
pathways_genes2 <- as.list(cgp_pathways)
paths <- as.data.frame(elementMetadata(cgp_pathways)[1])
names(pathways_genes2) <- paths$name
p <- pathways_genes2[grepl('BREAST',names(pathways_genes2))]
p <- append(p, pathways_genes2[grepl('ESR',names(pathways_genes2))])
pathfrac <- as.data.frame(names(pathways_genes2))
fraction <- c()
for (each in 1:length(names(pathways_genes2))){
  print(names(pathways_genes2)[each])
  print(pathways_genes2[[each]])
  f <- sum(pathways_genes2[[each]] %in% transbar$geneid)
  fraction <- append(fraction, f/length(pathways_genes2[[each]]))
}
pathfrac <- cbind(pathfrac, fraction)


##############################
# regression
regrbar <- read.csv('regr_GFPBC_libB_37798_feature_importance.csv')

geneIDs <- ensembldb::select(EnsDb.Hsapiens.v79, keys= regrbar$feature, keytype = "SYMBOL", columns = c("SYMBOL","GENEID"))
geneIDs <- geneIDs[!duplicated(geneIDs$SYMBOL),]
regrbar$geneid <- ifelse(regrbar$feature %in% geneIDs$SYMBOL, geneIDs$GENEID, regrbar$feature)

regrranks <- regrbar$fc_mae
names(regrranks) <- regrbar$geneid

# hallmark pathways
hallmark_pathways <- getFeatureSetCollection('GMTY42:human/H.gmt.bz2@REVISION-1')
pathways_genes <- as.list(hallmark_pathways)
paths <- as.data.frame(elementMetadata(hallmark_pathways)[1])
names(pathways_genes) <- paths$name
regr_fgsea_hall <- fgsea(pathways_genes, regrranks, minSize = 10, maxSize = 500)

##############################
#
# HEATMAP MAE 1
get_featmat <- function(keyword){
  files <- list.files(path=".", pattern="*.csv", full.names=TRUE, recursive=FALSE)
  bardf <- data.frame()
  for (f in files){
    if (grepl(keyword, f, fixed=TRUE)){
      bar <- read.csv(f)
      if (dim(bardf)[1] == 0){
        bardf <- as.data.frame(bar$fc_mae)
        row.names(bardf) <- bar$feature
        colnames(bardf)[dim(bardf)[2]] <- f
        print(bardf)
      }
      else{
        bar2 <- as.data.frame(bar[,-1:-2])
        row.names(bar2) <- bar$feature
        colnames(bar2) <- f
        bardf2 <- merge(bardf, bar2, by='row.names')
        remove(bardf)
        bardf <- as.data.frame(bardf2)
        row.names(bardf) <- bardf2$Row.names
        bardf <- bardf[,-1]
      }
    }
  }
  bardf
}

empirical_stats <- function(thedata){
  es <- c()
  for (v in thedata){
    if (v > 0){
      below <- sum(thedata <= -v)
      if (below == 0){
        es <- append(es, 0)
      }
      else{
        above <- sum(thedata >= v)
        if ((below/above) > 1){
          es <- append(es, 1)
        }
        else{
          es <- append(es, below/above)
        }
      }
    }
    else{
      es <- append(es, 1)
    }
  }
  es
}

statmat <- function(alldata){
  for (c in 1:length(colnames(alldata))){
    dat <- alldata[,c]
    s <- empirical_stats(dat)
    if (c == 1){
      estatmat <- as.data.frame(s)
    }
    else{
      estatmat <- cbind(estatmat, s)
    }
  }
  row.names(estatmat) <- row.names(alldata)
  colnames(estatmat) <- colnames(alldata)
  estatmat
}

getgenes <- function(df1){
  ids <- row.names(df1)
  annots <- select(org.Hs.eg.db, keys=ids, 
                   columns="SYMBOL", keytype="ENSEMBL")
  annots <- annots[!duplicated(annots$ENSEMBL), ]
  row.names(annots) <- annots$ENSEMBL
  row.names(df1) <- make.names(annots$SYMBOL, unique=TRUE)
  return(df1)
}


bardf1 <- get_featmat('transformer_1')
bardfstats1 <- statmat(bardf1)
bardf2 <- data.frame(lapply(bardf1, function(x) scale(x, center = FALSE, scale = max(x, na.rm = TRUE)/1)))
# sometimes needed
bardf1 <- getgenes(bardf1)
row.names(bardf2) <- row.names(bardf1)



# df[rowSums(df > value thresh) > in # of cols, ]
df <- bardf2[rowSums(bardf2 > 0.3) >= 3, ]


bars <- c('10678', '20182', '31581', '37798', '40758', '83026', '90850', '92196', '94397')
bars <- c('Depl', 'Surv')
ComplexHeatmap::Heatmap(df, name = 'Scaled Impact', 
                        row_names_gp = grid::gpar(fontsize = 8), 
                        column_names_gp = grid::gpar(fontsize = 8), 
                        column_labels = bars)

### 0.75 MAE ###
bardf0.75 <- get_featmat('transformer_0.75')
bardfstats0.75 <- statmat(bardf0.75)
bardf2 <- data.frame(lapply(bardf0.75, function(x) scale(x, center = FALSE, scale = max(x, na.rm = TRUE)/1)))
# sometimes needed
bardf0.75 <- getgenes(bardf0.75)
row.names(bardf2) <- row.names(bardf0.75)

# df[rowSums(df > value thresh) > in # of cols, ]
df <- bardf2[rowSums(bardf2 > 0.2) >= 2, ]


#bars <- c('10678', '20182', '31581', '37798', '40758', '83026', '90850', '92196', '94397')
ComplexHeatmap::Heatmap(df, name = 'Scaled Impact', 
                        row_names_gp = grid::gpar(fontsize = 4), 
                        column_names_gp = grid::gpar(fontsize = 4), 
                        column_labels = bars)



bardf0.5 <- get_featmat('transformer_0.5')
bardfstats0.5 <- statmat(bardf0.5)
bardf2 <- data.frame(lapply(bardf0.5, function(x) scale(x, center = FALSE, scale = max(x, na.rm = TRUE)/1)))
# sometimes needed
bardf0.5 <- getgenes(bardf0.5)
row.names(bardf2) <- row.names(bardf0.5)

# df[rowSums(df > value thresh) > in # of cols, ]
df <- bardf2[rowSums(bardf2 > 0.2) >= 2, ]


#bars <- c('10678', '20182', '31581', '37798', '40758', '83026', '90850', '92196', '94397')
ComplexHeatmap::Heatmap(df, name = 'Scaled Impact', 
                        row_names_gp = grid::gpar(fontsize = 4), 
                        column_names_gp = grid::gpar(fontsize = 4), 
                        column_labels = bars)

bardf0.25 <- get_featmat('transformer_0.25')
bardfstats0.25 <- statmat(bardf0.25)
bardf2 <- data.frame(lapply(bardf0.25, function(x) scale(x, center = FALSE, scale = max(x, na.rm = TRUE)/1)))
# sometimes needed
bardf0.25 <- getgenes(bardf0.25)
row.names(bardf2) <- row.names(bardf0.25)

# df[rowSums(df > value thresh) > in # of cols, ]
df <- bardf2[rowSums(bardf2 > 0.2) >= 2, ]


#bars <- c('10678', '20182', '31581', '37798', '40758', '83026', '90850', '92196', '94397')
ComplexHeatmap::Heatmap(df, name = 'Scaled Impact', 
                        row_names_gp = grid::gpar(fontsize = 4), 
                        column_names_gp = grid::gpar(fontsize = 4), 
                        column_labels = bars)


pathsig <- function(genenames, ranks, path){
  names(ranks) <- genenames
  pathways_genes <- as.list(path)
  paths <- as.data.frame(elementMetadata(path)[1])
  names(pathways_genes) <- paths$name
  pathscore <- fgsea(pathways_genes, ranks, minSize = 10, maxSize = 500)
  pathscore
}

sig <- list()
val <- 0.2
getgeneon <- FALSE
for (c in 1:dim(bardf1)[2]){
  siggenes <- c()
  x1 <- ifelse(bardfstats1[,c] <= val, row.names(bardfstats1), 'no')
  #test <- bardfstats1[bardfstats1 <= val, ]
  
  colgenes <- x1[x1 != 'no']
  #pathx1 <- pathsig(colgenes, bardfstats0.25[row.names(bardfstats0.25) %in% colgenes, c])
  
  
  siggenes <- append(siggenes, colgenes)
  x0.75 <- ifelse(bardfstats0.75[,c] <= val, row.names(bardfstats0.75), 'no')
  colgenes <- x0.75[x0.75 != 'no']
  siggenes <- append(siggenes, colgenes)
  x0.5 <- ifelse(bardfstats0.5[,c] <= val, row.names(bardfstats0.5), 'no')
  colgenes <- x0.5[x0.5 != 'no']
  siggenes <- append(siggenes, colgenes)
  x0.25 <- ifelse(bardfstats0.25[,c] <= val, row.names(bardfstats0.25), 'no')
  colgenes <- x0.25[x0.25 != 'no']
  siggenes <- append(siggenes, colgenes)
  siggenes <- unique(siggenes)
  print(siggenes)
  s <- strsplit(colnames(bardf1)[c], 'transformer_1')[[1]][2]
  s2 <- strsplit(s, '_feature')[[1]][1]
  print(s2)
  
  sigdf <- bardfstats0.25[row.names(bardfstats0.25) %in% siggenes, c]
  sigdf <- cbind(sigdf, bardfstats0.5[row.names(bardfstats0.5) %in% siggenes, c])
  sigdf <- cbind(sigdf, bardfstats0.75[row.names(bardfstats0.75) %in% siggenes, c])
  sigdf <- cbind(sigdf, bardfstats1[row.names(bardfstats1) %in% siggenes, c])
  print(sigdf)
  if (getgeneon == TRUE){
    ids <- siggenes
    annots <- select(org.Hs.eg.db, keys=ids, 
                     columns="SYMBOL", keytype="ENSEMBL")
    annots <- annots[!duplicated(annots$ENSEMBL), ]
    row.names(annots) <- annots$ENSEMBL
    siggenes <- make.names(annots$SYMBOL, unique=TRUE)
  }
  
  heat <- ComplexHeatmap::Heatmap(sigdf, name = 'Empirical Significance', 
                          row_names_gp = grid::gpar(fontsize = 6), 
                          column_names_gp = grid::gpar(fontsize = 10), 
                          row_labels = siggenes,
                          column_labels = c('0.25', '0.5', '0.75', '1'), 
                          cluster_columns = FALSE, column_title = s2)
  print(heat)
  
  for (c in 1:dim(sigdf)[2]){
    print(s2)
    print(sigdf[,c])
  }
}
  
