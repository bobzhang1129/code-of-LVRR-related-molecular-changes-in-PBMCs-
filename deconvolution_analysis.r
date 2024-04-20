## ----global_options, include=FALSE----------------------------------------------------------------------------------
library(knitr)
knitr::opts_chunk$set(dpi = 600, echo= TRUE, warning=FALSE, message=FALSE, fig.align = 'center',
                      fig.show=TRUE, fig.keep = 'all', fig.height= 6, fig.width=7)
    library(Seurat)
    library(dplyr)
    library(cowplot)
    library(ggplot2)
    library(pheatmap)
    library(rafalib)
    library(scPred)
    library(openxlsx)

library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(DescTools)
library(ggplot2)
library (VennDiagram)
library(reshape2)
library(gdata)
library(pheatmap)
library(PCAtools)

library(ReactomePA)
library(msigdbr)

library(enrichplot)
  library(Seurat)
  library(venn)
  library(dplyr)
  library(cowplot)
  library(ggplot2)
  library(pheatmap)
  library(enrichR)
  library(rafalib)



gwb <- createWorkbook()
###gsea

alldata <- readRDS("E:\\RNAseq\\scRNA\\pbmc_data\\pbmc_merge\\data\\results\\seurat_qc_dr_int_cl_dge_cell_type_remove.rds")
alldata
head(alldata@meta.data)
table(alldata@meta.data["orig.ident"])
table(alldata$orig.ident)
as.data.frame(alldata@assays$RNA@counts[1:10, 1:2])
print(names(alldata@reductions))
#str(alldata)
DefaultAssay(alldata)
###

sel.clust = "pathway"
alldata <- SetIdent(alldata, value = sel.clust)
table(alldata@active.ident)


markers_genes <- FindAllMarkers(alldata,
                               logfc.threshold = 0.2,
                               test.use = "wilcox",
                             min.pct = 0.1,
 #                              min.diff.pct = 0.2,
                               only.pos = TRUE,
                               max.cells.per.ident = 30,
#                               assay = "RNA")
                               assay = "CCA")


#markers_genes <- FindAllMarkers(alldata,
#                              logfc.threshold = 0.2,
#                               test.use = "LR",
 #                              latent.vars = 'Anatomical_Site',
 #                              min.pct = 0.1,
 #                              min.diff.pct = 0.2,
#                               only.pos = TRUE,
#                               max.cells.per.ident = 50,
#                               assay = "RNA")




#markers_genes <- FindAllMarkers(alldata,
#                              logfc.threshold = 0.2,
 ##                              test.use = "LR",
 ##                              latent.vars = 'orig.ident',
 #                              min.pct = 0.1,
 #                              min.diff.pct = 0.2,
 #                              only.pos = TRUE,
 #                              max.cells.per.ident = 50,
 #                              assay = "RNA")



markers_genes %>% group_by(cluster)  %>% top_n(-25, p_val_adj) -> top25
top25

#We can now select the top 25 up regulated genes for plotting.
mypar(2,5,mar=c(4,6,3,1))
for(i in unique(top25$cluster)){
  barplot( sort( setNames(top25$avg_log2FC, top25$gene) [top25$cluster == i], F),
           horiz = T,las=1 ,main=paste0(i," vs. rest"),border = "white", yaxs="i" )
  abline(v=c(0,0.25),lty=c(1,2))
}




####We can visualize them as a heatmap. Here we are selecting the top 5.
markers_genes %>% group_by(cluster)  %>% top_n(-5, p_val_adj) -> top5

# create a scale.data slot for the selected genes
alldata <- ScaleData(alldata, features = as.character(unique(top5$gene)), assay = "RNA")
DoHeatmap(alldata, features = as.character(unique(top5$gene)),group.by = sel.clust, assay = "RNA")


###
#Another way is by representing the overal group expression and detection rates in a dot-plot.

DotPlot(alldata, features = rev(as.character(unique(top5$gene))),group.by = sel.clust,assay = "RNA")+coord_flip()

# take top 3 genes per cluster/
top5 %>% group_by(cluster)  %>% top_n(-3, p_val) -> top3


# set pt.size to zero if you do not want all the points to hide the violin shapes, or to a small value like 0.1
VlnPlot(alldata, features = as.character(unique(top3$gene)), ncol = 5, group.by = sel.clust, assay = "RNA", pt.size = 0)

markers_genes <- markers_genes[!duplicated(markers_genes[c("cluster","gene")]),]

addWorksheet(gwb,sheetName = 'single_cell_mark')
writeData(gwb,sheet = 'single_cell_mark',markers_genes)

###bulid DB
DBi <- 50

DB <- markers_genes %>% group_by(cluster)  %>% top_n(-DBi, p_val) %>% filter(p_val < 0.01) %>% 
      dplyr::select(cluster,gene) %>% data.frame
#DB <- split(as.character(DB$gene), DB$cluster)

table(DB$cluster)
#lapply(DB,length)

### GENEs
genes <- read.xlsx("../GSEA_EA/transcriptomics_report.xlsx",sheet="all_genes") %>% 
         .[!duplicated(.["SYMBOL"]),] %>%filter( !is.na(SYMBOL))
#DEG_annotation

addWorksheet(gwb,sheetName = 'bluk_genes')
writeData(gwb,sheet = 'bluk_genes',genes)

genes <- genes %>% arrange(desc(logFC)) %>% filter(PValue<0.05)  %>% filter(abs(logFC) > log2(1.5))
#genes <- genes  %>% arrange(desc(logFC)) %>% filter(PValue<0.01) 
dim(genes)

geneList <- genes$logFC
names(geneList) <- genes$SYMBOL

GS <- GSEA(geneList, TERM2GENE= DB, pvalueCutoff = 1, minGSSize = 5, verbose=FALSE) %>% data.frame
GS


addWorksheet(gwb,sheetName = 'cell_type_GSEA')
writeData(gwb,sheet = 'cell_type_GSEA',GS)




##granulator

library(granulator)
#sel.clust <- "pathway"
#alldata@assays$RNA@data
#alldata@assays$RNA@scale.data


##diff
sel.clust = "pathway"
alldata <- SetIdent(alldata, value = sel.clust)
table(alldata@active.ident)


markers_genes <- FindAllMarkers(alldata,
                               logfc.threshold = 0.2,
                               test.use = "wilcox",
                             min.pct = 0.1,
 #                              min.diff.pct = 0.2,
                               only.pos = TRUE,
                               max.cells.per.ident = 30,
#                               assay = "RNA")
                               assay = "CCA")


#markers_genes <- markers_genes %>% filter(p_val < 0.05)




 sigMatrix <-   AverageExpression(
                alldata,
                assays = NULL,
                features = NULL,
                return.seurat = FALSE,
                group.by = sel.clust,
 #   add.ident = NULL,
#    slot = "data",
#    verbose = TRUE,
)

names(sigMatrix)

matrix_CCA <-  sigMatrix[[2]]


#matrix_CCA <- matrix_CCA[rownames(matrix_CCA) %in% rownames(markers_genes), ]

#matrix_CCA <- 2^(matrix_CCA)

max(matrix_CCA)
min(matrix_CCA)

matrix_CCA <- matrix_CCA - min(matrix_CCA)

max(matrix_CCA)
min(matrix_CCA)

dim(matrix_CCA)
#matrix_CCA <-  sigMatrix[[2]] %>% data.frame %>% as.matrix() 

#any(is.na(matrix_CCA))

#colnames(matrix_CCA) <- gsub("",".",colnames(matrix_CCA))

granulator::plot_similarity(sigMatrix=matrix_CCA)


logCPM <- read.table("../logCPM.xls",header = T,row.names=1)
counts <- 2^(logCPM) %>% as.matrix()

### GENEs
genes <- read.xlsx("../GSEA_EA/transcriptomics_report.xlsx",sheet="all_genes") %>% 
         .[!duplicated(.["SYMBOL"]),] %>%filter( !is.na(SYMBOL)) %>% 
         .[!duplicated(.["ENSEMBL"]),]


#genes <- genes %>% arrange(desc(logFC)) %>% filter(PValue<0.05)  %>% filter(abs(logFC) > log2(1.5))
genes <- genes %>% arrange(desc(logFC)) %>% filter(PValue<0.5)
dim(genes)
counts <- counts[genes$ENSEMBL,]
rownames(counts) <- genes$SYMBOL

max(counts)
min(counts)
dim(counts)
decon <- deconvolute(m = counts, sigMatrix = matrix_CCA)

plot_proportions(deconvoluted = decon, method = 'svr')


decon$proportions[["svr_sig1"]] %>% rowSums

plot_proportions(deconvoluted = decon, method = 'dtangle')
decon$proportions[[1]] %>% rowSums
##      B.Memory B.Naive Basophils.LD MAIT Monocytes.C
## 36TS     2.73    3.18         1.82 0.00       25.91
## 453W     2.27    4.09         0.45 0.45       24.09
## 4DUY     4.09    5.00         1.36 1.36       15.91
## 684C     1.36   11.82         1.36 1.82       21.82
## 925L     2.27   17.73         2.27 0.91       22.73
#We can plot the estimated cell type proportions with the function plot_proportions(). Notice that while the sum of cell types proportions cannot exceed 100%, for some methods part of the bulk RNA-seq signal remains unassigned.

# plot cell type proportions for svr model on ABIS_S0 reference profile
#plot_proportions(deconvoluted = decon, method = 'svr', signature = 'ABIS_S0')


#To plot all estimated cell type proportions we use the function plot_deconvolute(), which allows to compare results across deconvolution methods and cell types. The option scale indicates whether cell type proportions should be transformed into standard scores. Scaling is useful to directly compare deconvolution output, as the absolute percentages may vary considerably across methods.

# plot cell type proportions
plot_deconvolute(deconvoluted = decon, scale = TRUE, labels = FALSE)


decon_t <- decon$proportions[[1]]






















#################test
## ----global_options, include=FALSE----------------------------------------------------------------------------------
library(knitr)
knitr::opts_chunk$set(dpi = 600, echo= TRUE, warning=FALSE, message=FALSE, fig.align = 'center',
                      fig.show=TRUE, fig.keep = 'all', fig.height= 6, fig.width=7)
    library(Seurat)
    library(dplyr)
    library(cowplot)
    library(ggplot2)
    library(pheatmap)
    library(rafalib)
    library(scPred)
    library(openxlsx)

library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(DescTools)
library(ggplot2)
library (VennDiagram)
library(reshape2)
library(gdata)
library(pheatmap)
library(PCAtools)

library(ReactomePA)
library(msigdbr)

library(enrichplot)
  library(Seurat)
  library(venn)
  library(dplyr)
  library(cowplot)
  library(ggplot2)
  library(pheatmap)
  library(enrichR)
  library(rafalib)

library(MuSiC)
library(Biobase)
library(MuSiC2)
library(SingleCellExperiment)
library(scater)

setwd("E:\\RNAseq\\heart_failure\\mRNA\\celltype")


######MuSiC Deconvolution
logCPM <- read.table("../logCPM.xls",header = T,row.names=1)
counts <- 2^(logCPM)

### GENEs
genes <- read.xlsx("../GSEA_EA/transcriptomics_report.xlsx",sheet="all_genes") %>% 
         .[!duplicated(.["SYMBOL"]),] %>%filter( !is.na(SYMBOL)) %>% 
         .[!duplicated(.["ENSEMBL"]),]
#DEG_annotation

counts <- counts[genes$ENSEMBL,]
rownames(counts) <- genes$SYMBOL


###single cell
alldata <- readRDS("E:\\RNAseq\\scRNA\\pbmc_data\\pbmc_merge\\data\\results\\seurat_qc_dr_int_cl_dge_cell_type_remove.rds")
alldata
head(alldata@meta.data)
table(alldata@meta.data["orig.ident"])
table(alldata$orig.ident)
as.data.frame(alldata@assays$RNA@counts[1:10, 1:2])
print(names(alldata@reductions))
#str(alldata)
DefaultAssay(alldata)
###

sel.clust = "pathway"

table(alldata@meta.data[,sel.clust])

alldata <- SetIdent(alldata, value = alldata)
table(alldata@active.ident)

seger.sce <- as.SingleCellExperiment(alldata)


p1 <- plotExpression(seger.sce, features = "MS4A1", x = sel.clust) + theme(axis.text.x = element_text(angle = 45,
    hjust = 1))
p2 <- plotPCA(seger.sce, colour_by = sel.clust)
p1 + p2




prop_music=music_prop(bulk.mtx = counts, sc.sce = seger.sce, 
                     clusters = sel.clust, 
                     samples = 'orig.ident'
                     )


$Est.prop.weighted






####ref
clini <- read.table("../clini.xls",header = T,row.names=1)

bulk.control <- clini[clini$Symtomatic == "AS",] %>% rownames(.) 
bulk.control.mtx <-  counts[colnames(counts) %in% bulk.control]
bulk.control.mtx[1:5,1:5]

bulk.case <- clini[clini$Symtomatic == "S",] %>% rownames(.) 
bulk.case.mtx <-  counts[colnames(counts) %in% bulk.case]
bulk.case.mtx[1:5,1:5]



##
alldata <- readRDS("E:/RNAseq/scRNA/pbmc_data/pbmc10k/data/results/pbmc10k_qc_dr_int_cl_celltype.rds")
alldata
head(alldata@meta.data)
table(alldata@meta.data["orig.ident"])
table(alldata$orig.ident)
as.data.frame(alldata@assays$RNA@counts[1:10, 1:2])
print(names(alldata@reductions))
#str(alldata)
DefaultAssay(alldata)
###

sel.clust = "predicted.id" ####cell type
alldata <- SetIdent(alldata, value = sel.clust)
table(alldata@active.ident)

select_ct <- unique(GS$ID)
select_ct

alldata$Sample_ID <- sample (c(1,2,3), size=nrow(alldata@meta.data), replace=T)

#seger.sce <- SingleCellExperiment(list(counts=as.data.frame(alldata@assays$RNA@counts)),
#    colData=DataFrame(sampleID=alldata@meta.data[,"orig.ident"], cellType=alldata@meta.data[,sel.clust]),
#    rowData=DataFrame(SYBMOL= rownames(as.data.frame(alldata@assays$RNA@counts))),
#    metadata=list(study="PBMC10k")
#)

seger.sce <- as.SingleCellExperiment(alldata)


p1 <- plotExpression(seger.sce, features = "MS4A1", x = "ident") + theme(axis.text.x = element_text(angle = 45,
    hjust = 1))
p2 <- plotPCA(seger.sce, colour_by = "ident")
p1 + p2


head(colData(seger.sce))
assay(seger.sce)[1:5,1:5]

# music2 deconvolution
set.seed(1234)
est = music2_prop_t_statistics(bulk.control.mtx = bulk.control.mtx, 
                              bulk.case.mtx = bulk.case.mtx, sc.sce = seger.sce, 
                              clusters = "ident", 
                               samples = "Sample_ID", 
                              select.ct = select_ct, 
                              n_resample=20, sample_prop=0.5,cutoff_c=0.05,cutoff_r=0.01)



#est <- music_prop(bulk.mtx = bulk.control.mtx, sc.eset = seger.sce, clusters = "ident",
 #                              samples = 'orig.ident', select.ct = select_ct, verbose = F)




est.prop = est$Est.prop


# plot estimated cell type proportions
prop_all = cbind('proportion'=c(est.prop), 'sampleID'=rep(rownames(est.prop),times=ncol(est.prop)), 'celltype'=rep(colnames(est.prop), each=nrow(est.prop)))
prop_all = as.data.frame(prop_all)
prop_all$proportion = as.numeric(as.character(prop_all$proportion))
prop_all$group = ifelse(prop_all$sampleID %in% seq(from=1, to=100, by=1), 'Healthy', 'T2D')
cols <-c("alpha" = "cadetblue2", "beta" = "lightsalmon1", "delta" = "palegreen2", "ductal" = "goldenrod1",
          "gamma"="steelblue3", "acinar" = "plum2")
ggplot(prop_all, aes(x=celltype, y=proportion, color=celltype)) + xlab('')+
  geom_jitter(width=0.25,alpha=0.8)+ylab('Cell Type Proportions')+theme_bw()+
  stat_summary(fun = median,
               geom = "crossbar", width = 0.5,size=0.5,color='gray36')+
  facet_grid(.~group)+
  theme(plot.title = element_text(hjust = 0.5, size=12),
        axis.text.x = element_text(size=12,angle = 45,hjust=1),
        axis.text.y = element_text(size=12),
        axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12),
        axis.line = element_line(colour = "black"),
        strip.text.x = element_text(size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = 'none')+
  scale_color_manual(values=cols)







#saveWorkbook(gwb, "transcriptomics_report.xlsx", overwrite = TRUE)