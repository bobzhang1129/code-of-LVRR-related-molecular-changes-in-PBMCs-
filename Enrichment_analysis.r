## ----global_options, include=FALSE----------------------------------------------------------------------------------
library(knitr)
knitr::opts_chunk$set(dpi = 300, echo= TRUE, warning=FALSE, message=FALSE, fig.align = 'center',
                      fig.show=TRUE, fig.keep = 'all', fig.height= 6, fig.width=7)
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)v
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
library(ReactomePA)
library(msigdbr)
library(dplyr)
library(curl)
library(openxlsx)
library(enrichplot)


#GSEA
bp <- gseGO(geneList     = geneList,
              OrgDb        = org.Hs.eg.db,
              ont          = "BP",
              minGSSize    = 10,
              pvalueCutoff = 1,
#	      pAdjustMethod = "none",
              verbose      = FALSE)

bp <- setReadable(bp, 'org.Hs.eg.db', keyType="ENTREZID")



cc <- gseGO(geneList     = geneList,
              OrgDb        = org.Hs.eg.db,
              ont          = "BP",
              minGSSize    = 10,
              pvalueCutoff = 1,
#	      pAdjustMethod = "none",
              verbose      = FALSE)

cc <- setReadable(cc, 'org.Hs.eg.db', keyType="ENTREZID")



mf <- gseGO(geneList     = geneList,
              OrgDb        = org.Hs.eg.db,
              ont          = "BP",
              minGSSize    = 10,
              pvalueCutoff = 1,
#	      pAdjustMethod = "none",
              verbose      = FALSE)

mf <- setReadable(mf, 'org.Hs.eg.db', keyType="ENTREZID")



#########EA


bp <- enrichGO(gene     = names(geneList),
              OrgDb        = org.Hs.eg.db,
              ont          = "BP",
              minGSSize    = 10,
              maxGSSize = 3000,
              pvalueCutoff = p_sel,
              qvalueCutoff  = 1,
#	      pAdjustMethod = "none",
              )

bp <- setReadable(bp, 'org.Hs.eg.db', keyType="ENTREZID")


cc <- enrichGO(gene     = names(geneList),
              OrgDb        = org.Hs.eg.db,
              ont          = "CC",
              minGSSize    = 10,
              maxGSSize = 3000,
              pvalueCutoff = 1,
              qvalueCutoff  = 1,
#	      pAdjustMethod = "none",
              )

cc <- setReadable(cc, 'org.Hs.eg.db', keyType="ENTREZID")

mf <- enrichGO(gene     = names(geneList),
              OrgDb        = org.Hs.eg.db,
              ont          = "MF",
              minGSSize    = 10,
              maxGSSize = 3000,
              pvalueCutoff = 1,
              qvalueCutoff  = 1,
#	      pAdjustMethod = "none",
              )

mf <- setReadable(mf, 'org.Hs.eg.db', keyType="ENTREZID")

