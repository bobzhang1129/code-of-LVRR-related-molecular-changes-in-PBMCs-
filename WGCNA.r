## ----global_options, include=FALSE----------------------------------------------------------------------------------
library(knitr)
knitr::opts_chunk$set(dpi = 600, echo= TRUE, warning=FALSE, message=FALSE, fig.align = 'center',
                      fig.show=TRUE, fig.keep = 'all', fig.height= 6, fig.width=7)
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
library(gtools)
library(pathview)
library(ggplot2)
library(stringr)
library(GOplot)
library(dplyr)
library(pathview)
library(KEGGREST)
library(dplyr)
library(DescTools)
library(ggplot2)
library(showtext)
library(psych)
library(corrplot)
library(stringr)
library (VennDiagram)
library(reshape2)
library(gdata)
library(pheatmap)
library(PCAtools)
library(kmodR)
library(ropls)
library(limma)
library(statmod)
library(openxlsx)
library(ReactomePA)
library(msigdbr)
library(dplyr)
library(curl)
library(enrichplot)
library(plotrix)
library(WGCNA)


setwd(workingDir)
options(stringsAsFactors = FALSE) 


gsg = goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK

if (!gsg$allOK)
{
  
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
 
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

sampleTree = hclust(dist(datExpr0), method = "average")
sizeGrWindow(12,9) 
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)


cf <-85

abline(h = cf, col = "red") 
clust = cutreeStatic(sampleTree, cutHeight = cf, minSize = 10)

table(clust)   
keepSamples = (clust==1) 
datExpr = datExpr0[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)


datTraits <- metadata[-1:-2]
datTraits$group <- ifelse(datTraits$group == "LVRR",1,0)

sampleTree2 = hclust(dist(datExpr), method = "average")
traitColors = numbers2colors(datTraits, signed = FALSE) #用颜色代表关联度
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap")


save(datExpr, datTraits, file = "WGCNA-dataInput.RData")


getwd();

setwd(workingDir);
library(WGCNA)
library(openxlsx)
options(stringsAsFactors = FALSE);
enableWGCNAThreads() 
lnames = load(file = "WGCNA-dataInput.RData");
lnames

powers = c(c(1:10), seq(from = 12, to=24, by=2))
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;

plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
     abline(h=0.90,col="red")  #查看位于0.9以上的点，可以改变高度值


plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

sft$powerEstimate


power <- sft$powerEstimate

type <- "signed"

if (is.na(power)){ 
  power = ifelse(nSamples<20, ifelse(type == "unsigned", 9, 18), 
                 ifelse(nSamples<30, ifelse(type == "unsigned", 8, 16), 
                        ifelse(nSamples<40, ifelse(type == "unsigned", 7, 14), 
                               ifelse(type == "unsigned", 6, 12))        
                 ) 
  ) 
}
power <- 20


softPower = power
adjacency = adjacency(datExpr, power = softPower)

TOM = TOMsimilarity(adjacency)

dissTOM = 1 - TOM

geneTree = hclust(as.dist(dissTOM), method = "average")

sizeGrWindow(12,9)
plot(geneTree, xlab = "", sub = "", main = "Gene clustering on TOM-based dissimilarity",
    labels = FALSE, hang = 0.04)



minModuleSize = 30

dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = 2, pamRespectsDendro = FALSE,
    minClusterSize = minModuleSize)

table(dynamicMods)


dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)

sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", dendroLabels = FALSE,
    hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "Gene dendrogram and module colors")



MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes

MEDiss = 1 - cor(MEs)

METree = hclust(as.dist(MEDiss), method = "average")

sizeGrWindow(7,6)
plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")
MEDissThres = 0.3

abline(h = MEDissThres, col = "red")



merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)

mergedColors = merge$colors

mergedMEs = merge$newMEs
sizeGrWindow(12,9)

plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), c("Dynamic Tree Cut",
    "Merged dynamic"), dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)






moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder) - 1
MEs = mergedMEs
# Save module colors and labels for use in subsequent parts

table(moduleColors)
#library(dplyr)
#mdf <- data.frame(genes = colnames(datExpr),modules = moduleColors) 
#       %>% filter (genes %in% spgene)
#mdf
#table(mdf$modules)
#table(moduleColors)

save(MEs, moduleLabels, moduleColors, geneTree,power, file = "WGCNA-02-networkConstruction-stepByStep.RData")
