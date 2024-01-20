setwd('/gstore/project/hr_brca_heterogeneity/T47D_trace_Seq_v2')
sce <- gp.sa.core::readResult("traceseq_pathway_score-results/sce", ".", rooted=TRUE)

saveRDS(sce, '/gstore/project/hr_brca_heterogeneity/T47D_trace_Seq_v2/Nathan/nathan_sce.rds')
sce <- readRDS('/gstore/project/hr_brca_heterogeneity/T47D_trace_Seq_v2/Nathan/nathan_sce.rds')
