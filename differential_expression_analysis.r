## ----global_options, include=FALSE----------------------------------------------------------------------------------
library(knitr)
knitr::opts_chunk$set(dpi = 300, echo= TRUE, warning=FALSE, message=FALSE, fig.align = 'center',
                      fig.show=TRUE, fig.keep = 'all', fig.height= 6, fig.width=7)

library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
library(gtools)
library(pathview)
library(stringr)
library(GOplot)
library(ggplot2)
library(ggfortify)
library(cluster)
library(pheatmap)
library(randomcoloR)
library(edgeR)
library(DESeq2)
library(dendextend)
library(pathview)
library(gtools)
library(KEGGREST)
library (VennDiagram)
library(ggrepel)
library(dplyr)
library(openxlsx)
library(PCAtools) 
#library(EnsDb.Mmusculus.v79)
library(rtracklayer)
library(ggbreak)
library(ggforce)
library(plotrix)
library(ggrepel)
#library(RcmdrMisc)
#library(psych)

y <- DGEList(counts=counts,group=clini$group)
keep <- filterByExpr(y)

y <- y[keep,,keep.lib.sizes=FALSE]

y <- calcNormFactors(y)

y$samples

logCPM <- cpm(y,log=TRUE)
write.table(logCPM, file="logCPM.xls", sep="\t", col.names=NA)

#
dim(counts) ##before QC
dim(logCPM) ##after QC




#####
cols <- as.numeric(y$samples$group)
plotMDS(y,col=cols)

###
mm <- model.matrix(~ group,clini)
#mm <- model.matrix(~ group,clini)
mm

#cpm<- removeBatchEffect(cpm, batch=ex, design=mm)
#plotMDS(cpm, col=cols)

#logCPM<- removeBatchEffect(logCPM, batch=ex, design=mm)
#plotMDS(logCPM, col=cols)
#write.table(logCPM, file="logCPM_pca.xls", sep="\t", col.names=NA)

#par(mfrow=c(2,1))
#plotMDS(cpm, col=cols)
#plotMDS(logCPM, col=cols)



y <- estimateGLMCommonDisp(y,mm)
y <- estimateGLMTrendedDisp(y,mm)
### Loading required package: splines
y <- estimateGLMTagwiseDisp(y,mm)
plotBCV(y) ##红线在0.2~0.4之间，否则需要考虑实验样品之间的差异是否过大了。


f <- glmFit(y,mm)

lrt_ab <- glmLRT(f, coef= 2)

topTags(lrt_ab, 20)
deall<- topTags(lrt_ab, n=Inf)$table

########annotation by coding noncoding
#gtf <- importGTF(file="E:\\RNAseq\\pcsk6ko\\QCdiff\\Mus_musculus.GRCm39.104.gtf.gz")
#edb <- EnsDb.Mmusculus.v79
#edb
#organism(edb)
#supportedFilters(edb)
#Tx <- transcripts(edb, filter = GeneNameFilter("BCL2L11"))
#transcripts(edb, filter = ~ gene_name == "ZBTB16")
#tx <- data.frame(transcripts(edb))
#rownames(tx) <- tx$gene_id
#dim(tx[tx$gene_id %in% rownames(deall),])

#gtf <- rtracklayer::import("E:\\RNAseq\\pcsk6ko\\QCdiff\\Mus_musculus.GRCm39.104.gtf")
#gtf_df=data.frame(gtf)
#gtf_df <- gtf_df[!is.na(gtf_df$gene_biotype),]
#gtf_df <- gtf_df[!duplicated(gtf_df$gene_id),]
#deall$gene_id <- rownames(deall)
#deall <- left_join(deall,gtf_df,by="gene_id")
#rownames(deall)<- deall$gene_id

write.table(deall, file = "diff_all.xls", sep="\t", col.names=NA)
