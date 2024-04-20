###R4.2+

#remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
#if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
#BiocManager::install("celldex")
#BiocManager::install("SingleR")
#BiocManager::install("glmGamPoi")


##Removal of ambient RNA using SoupX

#https://cellgeni.github.io/notebooks/html/new-10kPBMC-SoupX.html
#####soupX remove free RNA out of cells (background mRNAs)
#In droplet based, single cell RNA-seq experiments, 
#there is always a certain amount of background mRNAs present in the 
#dilution that gets distributed into the droplets with cells and sequenced along with them. 
#install.packages("SoupX")

setwd("E:\\RNAseq\\scRNA\\pbmc_data\\pbmc_merge")

library(Seurat)
library(SoupX)
library(DropletUtils)
library(ggplot2)
library(DoubletFinder)
library(knitr)
library(hdf5r)



filt.matrix <- Seurat::Read10X_h5("Parent_NGSC3_DI_PBMC_filtered_feature_bc_matrix.h5",use.names = T)

raw.matrix  <- Seurat::Read10X_h5("Parent_NGSC3_DI_PBMC_raw_feature_bc_matrix.h5",use.names = T)


str(raw.matrix)

str(filt.matrix)

srat  <- CreateSeuratObject(counts = filt.matrix,project = "PBMC")
srat$type = "srat"
srat

soup.channel  <- SoupChannel(raw.matrix, filt.matrix)
soup.channel


srat    <- SCTransform(srat, verbose = F)
srat    <- RunPCA(srat, verbose = F)
srat    <- RunUMAP(srat, dims = 1:30, verbose = F)
srat    <- FindNeighbors(srat, dims = 1:30, verbose = F)
srat    <- FindClusters(srat, verbose = T)

meta    <- srat@meta.data
umap    <- srat@reductions$umap@cell.embeddings
soup.channel  <- setClusters(soup.channel, setNames(meta$seurat_clusters, rownames(meta)))
soup.channel  <- setDR(soup.channel, umap)
head(meta)

soup.channel  <- autoEstCont(soup.channel)

head(soup.channel$soupProfile[order(soup.channel$soupProfile$est, decreasing = T), ], n = 20)

adj.matrix  <- adjustCounts(soup.channel, roundToInt = T)

DropletUtils:::write10xCounts("soupX_pbmc10k_filt", adj.matrix)


setwd("E:\\RNAseq\\scRNA\\pbmc_data\\pbmc_merge")

library(Seurat)
library(SoupX)
library(DropletUtils)
library(ggplot2)
library(DoubletFinder)
library(knitr)
library(hdf5r)
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(tidyr)
library(SingleR)
library(celldex)
library(RColorBrewer)
library(SingleCellExperiment)


setwd("E:\\RNAseq\\scRNA\\pbmc_data\\pbmc_merge")

adj.matrix <- Read10X("soupX_pbmc10k_filt")

srat <- CreateSeuratObject(adj.matrix,project = "pbmc10k") 
srat$type <- "pbmc10k"

str(srat)



pbmc.data <- Read10X(data.dir = "pbmc3k_filtered_gene_bc_matrices")

pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc
pbmc$type <- "pbmc3k"

head(pbmc@meta.data)
as.data.frame(pbmc@assays$RNA@counts[1:10, 1:2])


##5k


filt.matrix <- read.table("HPA\\read_count.tsv",header=T,sep="\t",row.names=1)

filt.matrix[1:5,1:5]

str(filt.matrix)

#all$ENSEMBL <- as.character(all$ENSEMBL)

order <- rowMeans(filt.matrix) %>% sort(decreasing = T)

filt.matrix <- filt.matrix[names(order),]

id_tran <- bitr(rownames(filt.matrix), fromType = "ENSEMBL",
         toType = c("SYMBOL"),
         OrgDb = org.Hs.eg.db)

id_tran <- id_tran %>% .[!duplicated(.$ENSEMBL),] %>% .[!duplicated(.$SYMBOL),]

filt.matrix <- filt.matrix[id_tran$ENSEMBL,]
rownames(filt.matrix) <- id_tran$SYMBOL

filt.matrix[1:5,1:5]

#str(filt.matrix)


#pbmc.data <- Read10X(data.dir = "pbmc3k_filtered_gene_bc_matrices")
# Initialize the Seurat object with the raw (non-normalized data).
pbmc_5k <- CreateSeuratObject(counts = filt.matrix, project = "pbmc5k", min.cells = 3, min.features = 200)
pbmc_5k
pbmc_5k$type <- "pbmc5k"

head(pbmc_5k@meta.data)
as.data.frame(pbmc_5k@assays$RNA@counts[1:10, 1:2])


#srat <- merge(srat, c(pbmc,pbmc_5k), add.cell.ids = c("pbmc10k", "pbmc3k","pbmc5k"))

srat <- merge(pbmc_5k, c(pbmc,srat), add.cell.ids = c("pbmc5k", "pbmc3k","pbmc10k"))

meta <- srat@meta.data
dim(meta)

head(meta)

summary(meta$nCount_RNA)

summary(meta$nFeature_RNA)

as.data.frame(srat@assays$RNA@counts[1:10, 1:2])
head(srat@meta.data, 10)

#Let¡¯s add several more values useful in diagnostics of cell quality. 
#Michochondrial genes are useful indicators of cell state. For mouse datasets, 
#change pattern to ¡°Mt-¡±, or explicitly list gene IDs with the features = ¡­ option.

#srat[["percent.mt"]] <- PercentageFeatureSet(srat, pattern = "^MT-")
srat <- PercentageFeatureSet(srat, "^MT-", col.name = "percent_mito")

#Similarly, we can define ribosomal proteins (their names begin with RPS or RPL), which often take substantial fraction of reads:

#srat[["percent.rb"]] <- PercentageFeatureSet(srat, pattern = "^RP[SL]")
srat <- PercentageFeatureSet(srat, "^RP[SL]", col.name = "percent_ribo")


#And finally, with the same method we will calculate proportion hemoglobin genes, 
#which can give an indication of red blood cell contamination.

# Percentage hemoglobin genes - includes all genes starting with HB except HBP.
srat <- PercentageFeatureSet(srat, "^HB[^(P)]", col.name = "percent_hb")

srat <- PercentageFeatureSet(srat, "PECAM1|PF4", col.name = "percent_plat")

head(srat@meta.data, 10)


##Plot QC

feats <- c("nFeature_RNA", "nCount_RNA", "percent_mito", "percent_ribo", "percent_hb")
VlnPlot(srat, group.by = "orig.ident", features = feats, pt.size = 0.1, ncol = 3) + NoLegend()


#Let¡¯s plot some of the metadata features against each other and see how they correlate. 
#The number above each plot is a Pearson correlation coefficient.

FeatureScatter(srat, "nCount_RNA", "nFeature_RNA", group.by = "orig.ident", pt.size = 0.5)

FeatureScatter(srat, feature1 = "nCount_RNA", feature2 = "percent_mito")

#Filtering
#Detection-based filtering


#A standard approach is to filter cells with low amount of reads as well as genes that are present in at least a certain amount of cells. 
#Here we will only consider cells with at least 200 detected genes and genes need to be expressed in at least 3 cells. 
#Please note that those values are highly dependent on the library preparation method used.

selected_c <- WhichCells(srat, expression = nFeature_RNA > 200)
selected_f <- rownames(srat)[Matrix::rowSums(srat) > 3]

data.filt <- subset(srat, features = selected_f, cells = selected_c)
dim(data.filt)
dim(srat)


#Additionally, we can also see which genes contribute the most to such reads. 
#We can for instance plot the percentage of counts per gene.

# Compute the relative expression of each gene per cell Use sparse matrix
# operations, if your dataset is large, doing matrix devisions the regular way
# will take a very long time.
par(mar = c(4, 8, 2, 1))
C <- data.filt@assays$RNA@counts
C <- Matrix::t(Matrix::t(C)/Matrix::colSums(C)) * 100
most_expressed <- order(apply(C, 1, median), decreasing = T)[20:1]
boxplot(as.matrix(t(C[most_expressed, ])), cex = 0.1, las = 1, xlab = "% total count per cell",
    col = (scales::hue_pal())(20)[20:1], horizontal = TRUE)

par(mar=c(5.1, 4.1, 4.1, 2.1))

#As you can see, MALAT1 constitutes up to 30% of the UMIs from a single cell and 
#the other top genes are mitochondrial and ribosomal genes. 
#It is quite common that nuclear lincRNAs have correlation with quality and mitochondrial reads, 
#so high detection of MALAT1 may be a technical issue. Let us assemble some information about such genes, 
#which are important for quality control and downstream filtering.


#Mito/Ribo filtering

#We also have quite a lot of cells with high proportion of mitochondrial and low proportion ofribosomal reads. It could be wise to remove those cells, if we have enough cells left after filtering.
#Another option would be to either remove all mitochondrial reads from the dataset and hope that the remaining genes still have enough biological signal.
#A third option would be to just regress out the percent_mito variable during scaling. In this case we had as much as 99.7% mitochondrial reads in some of the cells, so it is quite unlikely that there is much cell type signature left in those.
#Looking at the plots, make reasonable decisions on where to draw the cutoff. In this case, the bulk of the cells are below 20% mitochondrial reads and that will be used as a cutoff. We will also remove cells with less than 5% ribosomal reads.

head(data.filt@meta.data, 10)

summary(data.filt@meta.data["percent_mito"])
selected_mito <- WhichCells(data.filt, expression = percent_mito < 20)
length(selected_mito)


summary(data.filt@meta.data["percent_ribo"])
selected_ribo <- WhichCells(data.filt, expression = percent_ribo > 5)
length(selected_ribo)


# and subset the object to only keep those cells
data.filt <- subset(data.filt, cells = selected_mito)
data.filt <- subset(data.filt, cells = selected_ribo)

dim(data.filt)

table(data.filt$orig.ident)


#Plot filtered QC
#Lets plot the same QC-stats another time.

feats <- c("nFeature_RNA", "nCount_RNA", "percent_mito", "percent_ribo", "percent_hb")

VlnPlot(data.filt, group.by = "orig.ident", features = feats, pt.size = 0.1, ncol = 3) + NoLegend()


##Filter genes
#As the level of expression of mitochondrial and MALAT1 genes are judged as mainly technical, 
#it can be wise to remove them from the dataset bofore any further analysis.


dim(data.filt)

# Filter MALAT1
data.filt <- data.filt[!grepl("MALAT1", rownames(data.filt)), ]

# Filter Mitocondrial
data.filt <- data.filt[!grepl("^MT-", rownames(data.filt)), ]

# Filter Ribossomal gene (optional if that is a problem on your data) data.filt
# <- data.filt[ ! grepl('^RP[SL]', rownames(data.filt)), ]

# Filter Hemoglobin gene (optional if that is a problem on your data)
#data.filt <- data.filt[!grepl("^HB[^(P)]", rownames(data.filt)), ]

dim(data.filt)


##Sample sex (skip) check the detail of tuital


###Calculate cell-cycle scores


# Before running CellCycleScoring the data need to be normalized and
# logtransformed.
data.filt = NormalizeData(data.filt)


data.filt <- CellCycleScoring(object = data.filt, g2m.features = cc.genes$g2m.genes,
    s.features = cc.genes$s.genes)


VlnPlot(data.filt, features = c("S.Score", "G2M.Score"), group.by = "orig.ident",
    ncol = 4, pt.size = 0.1)



####Predict doublets

# remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
suppressMessages(require(DoubletFinder))


data.filt = FindVariableFeatures(data.filt, verbose = F)
data.filt = ScaleData(data.filt, vars.to.regress = c("nFeature_RNA", "percent_mito"),
    verbose = F)
data.filt = RunPCA(data.filt, verbose = F, npcs = 20)
data.filt = RunUMAP(data.filt, dims = 1:10, verbose = F)


######
#scDblFinder (R) for  Doublet identification   

#To optimize the parameters, you can run the paramSweep function in the package.

###find pk
# load libraries
library(parallel) # detectCores()
library(DoubletFinder) # paramSweep_v3()
library(Seurat) # SplitObject()


set.seed(8)

# work in parallel
options(mc.cores = detectCores() - 1)

table(data.filt$orig.ident)

filt.split <- SplitObject(data.filt, split.by = "orig.ident") 

# loop through samples to find doublets
for (i in 1:length(filt.split)) {
#  i=1
  # print the sample we are on
  print(paste0("orig.ident ",i))
  
  # Pre-process seurat object with standard seurat workflow
  data.sample <- NormalizeData(filt.split[[i]])
  data.sample <- FindVariableFeatures(data.sample)
  data.sample <- ScaleData(data.sample)
  data.sample <- RunPCA(data.sample, nfeatures.print = 10)
  
# Find significant PCs
  stdv <- data.sample[["pca"]]@stdev
  sum.stdv <- sum(data.sample[["pca"]]@stdev)
  percent.stdv <- (stdv / sum.stdv) * 100
  cumulative <- cumsum(percent.stdv)
  co1 <- which(cumulative > 90 & percent.stdv < 5)[1]
  co2 <- sort(which((percent.stdv[1:length(percent.stdv) - 1] - 
                       percent.stdv[2:length(percent.stdv)]) > 0.1), 
              decreasing = T)[1] + 1
  min.pc <- min(co1, co2)
  min.pc


  # finish pre-processing
  data.sample <- RunUMAP(data.sample, dims = 1:min.pc)
  data.sample <- FindNeighbors(object = data.sample, dims = 1:min.pc)              
  data.sample <- FindClusters(object = data.sample, resolution = 0.1)
  
  # pK identification (no ground-truth)
  #sweep.list <- paramSweep_v3(data.sample, PCs = 1:min.pc, num.cores = detectCores() - 1) ###linux
  sweep.list <- paramSweep_v3(data.sample, PCs = 1:min.pc, num.cores = )
  sweep.stats <- summarizeSweep(sweep.list)
  bcmvn <- find.pK(sweep.stats)
  
  # Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
  bcmvn.max <- bcmvn[which.max(bcmvn$BCmetric),]
  optimal.pk <- bcmvn.max$pK
  optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]
  
  ## Homotypic doublet proportion estimate
  annotations <- data.sample@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations) 
  nExp.poi <- round(optimal.pk * nrow(data.sample@meta.data)) ## Assuming 7.5% doublet formation rate - tailor for your dataset
  nExp.poi.adj <- round(nExp.poi * (1 - homotypic.prop))
  
  # run DoubletFinder
  data.sample <- doubletFinder_v3(seu = data.sample, 
                                   PCs = 1:min.pc, 
                                   pK = optimal.pk,
                                   nExp = nExp.poi.adj)
  metadata <- data.sample@meta.data
  colnames(metadata)[15] <- "DF.name"
  data.sample@meta.data <- metadata 
  
  # subset and save
  data.singlets <- subset(data.sample, DF.name == "Singlet")
  filt.split[[i]] <- data.singlets
  remove(data.singlets)
}


# converge filt.split  ##3 samples
#data.filt <- Reduce(function(x,y) merge(x,y,add.cell.ids = c(x@project.name,y@project.name)) , filt.split)


data.singlets <- merge(x = filt.split[[1]],
                       y = c(filt.split[[2]],filt.split[[3]]),
                       project = "PBMC scRNAseq")
data.singlets


data.filt


data.filt <- data.singlets

table(data.filt$orig.ident)

# name of the DF prediction can change, so extract the correct column name.
#DF.name = colnames(data.filt@meta.data)[grepl("DF.name", colnames(data.filt@meta.data))]



#cowplot::plot_grid(ncol = 2, DimPlot(data.filt, group.by = "orig.ident") + NoAxes(),
#    DimPlot(data.filt, group.by = DF.name) + NoAxes())


#VlnPlot(data.filt, features = "nFeature_RNA", group.by = DF.name, pt.size = 0.1)

#data.filt = data.filt[, data.filt@meta.data[, DF.name] == "Singlet"]
#dim(data.filt)

dir.create("data/results", showWarnings = T,recursive = TRUE)

saveRDS(data.filt, "data/results/seurat_qc.rds")



######Dimensionality reduction


suppressPackageStartupMessages({
    library(Seurat)
    library(cowplot)
    library(ggplot2)
    library(scran)
})

setwd("E:\\RNAseq\\scRNA\\pbmc_data\\pbmc_merge")

alldata <- readRDS("data/results/seurat_qc.rds")

head(alldata@meta.data)
as.data.frame(alldata@assays$RNA@counts[1:10, 1:2])

table(alldata@meta.data$orig.ident)

#Feature selection
#Next, we first need to define which features/genes are important in our dataset to distinguish cell types. 
#For this purpose, we need to find genes that are highly variable across cells, 
#which in turn will also provide a good separation of the cell clusters.

suppressWarnings(suppressMessages(alldata <- FindVariableFeatures(alldata, selection.method = "vst",
    nfeatures = 2000, verbose = FALSE, assay = "RNA")))
top20 <- head(VariableFeatures(alldata), 20)

LabelPoints(plot = VariableFeaturePlot(alldata), points = top20, repel = TRUE)


#we can use regression to remove any unwanted sources of variation from the dataset, such as cell cycle, 
#sequencing depth, percent mitocondria. 
#This is achieved by doing a generalized linear regression using these parameters as covariates 

alldata <- ScaleData(alldata, vars.to.regress = c("percent_mito", "nFeature_RNA"),
    assay = "RNA")

#Performing PCA has many useful applications and interpretations, 
#which much depends on the data used. In the case of life sciences, 
#we want to segregate samples based on gene expression patterns in the data.

alldata <- RunPCA(alldata, npcs = 50, verbose = F)

plot_grid(nrow = 3, DimPlot(alldata, reduction = "pca", group.by = "orig.ident",
    dims = 1:2), DimPlot(alldata, reduction = "pca", group.by = "orig.ident", dims = 3:4),
    DimPlot(alldata, reduction = "pca", group.by = "orig.ident", dims = 5:6))

plot_grid(nrow = 2, DimPlot(alldata, reduction = "pca", group.by = "orig.ident",
    dims = 1:2), DimPlot(alldata, reduction = "pca", group.by = "orig.ident", dims = 3:4))



#To identify which genes (Seurat) or metadata paramters (Scater/Scran) contribute the most to each PC, 
#one can retreive the loading matrix information. Unfortunatelly this is not implemented in Scater/Scran, 
#so you will need to compute PCA using logcounts.

VizDimLoadings(alldata, dims = 1:3, reduction = "pca", ncol = 3, balanced = T)

ElbowPlot(alldata, reduction = "pca", ndims = 50)

#Based on this plot, we can see that the top 8 PCs retain a lot of information, 
#while other PCs contain pregressivelly less. 
#However, it is still advisable to use more PCs since they might contain informaktion 
#about rare cell types (such as platelets and DCs in this dataset)


###tSNE

alldata <- RunTSNE(alldata, reduction = "pca", dims = 1:30, perplexity = 30, max_iter = 1000,
    theta = 0.5, eta = 200, num_threads = 0)
# see ?Rtsne and ?RunTSNE for more info

plot_grid(ncol = 1, DimPlot(alldata, reduction = "tsne", group.by = "orig.ident"))



###UMAP

alldata <- RunUMAP(alldata, reduction = "pca", dims = 1:30, n.components = 2, n.neighbors = 30,
    n.epochs = 200, min.dist = 0.3, learning.rate = 1, spread = 1)
# see ?RunUMAP for more info


#Your turn

#We have now done Variable gene selection, 
#PCA and UMAP with the settings we chose. Test a few different ways of selecting variable genes, 
#number of PCs for UMAP and check how it influences your embedding.


# we can add in additional reductions, by defulat they are named 'pca', 'umap',
# 'tsne' etc. But we can specify alternative names with reduction.name

alldata <- RunUMAP(alldata, reduction.name = "UMAP10_on_PCA", reduction = "pca",
    dims = 1:30, n.components = 10, n.neighbors = 30, n.epochs = 200, min.dist = 0.3,
    learning.rate = 1, spread = 1)
# see ?RunUMAP for more info

#We can now plot the UMAP colored per dataset. Although less distinct as in the tSNE, 
#we still see quite an effect of the different batches in the data.

plot_grid(nrow = 3, DimPlot(alldata, reduction = "umap", group.by = "orig.ident") +
    ggplot2::ggtitle(label = "UMAP_on_PCA"), DimPlot(alldata, reduction = "UMAP10_on_PCA",
    group.by = "orig.ident", dims = 1:2) + ggplot2::ggtitle(label = "UMAP10_on_PCA"),
    DimPlot(alldata, reduction = "UMAP10_on_PCA", group.by = "orig.ident", dims = 3:4) +
        ggplot2::ggtitle(label = "UMAP10_on_PCA"))



#We can now plot PCA, UMAP and tSNE side by side for comparison. 
#Here, we can conclude that our dataset contains a batch effect that needs to be corrected before 
#proceeding to clustering and differential gene expression analysis.

plot_grid(nrow = 3, DimPlot(alldata, reduction = "pca", group.by = "orig.ident"),
    DimPlot(alldata, reduction = "tsne", group.by = "orig.ident"), DimPlot(alldata,
        reduction = "umap", group.by = "orig.ident"))



###Using ScaledData and graphs for DR

#We will show from now on only UMAP, but the same applies for tSNE.

#####Using ScaledData for UMAP


#To run tSNE or UMAP on the scaled data, 
#one firts needs to select the number of variables to use. 



alldata <- RunUMAP(alldata, reduction.name = "UMAP_on_ScaleData", features = alldata@assays$RNA@var.features,
    assay = "RNA", n.components = 2, n.neighbors = 50, n.epochs = 200, min.dist = 0.3,
    learning.rate = 1, spread = 1)


#Build Graph
#To run tSNE or UMAP on the a graph, we first need to build a graph from the data. 
#In fact, both tSNE and UMAP first build a graph from the data using a specified distance metrix 
#and then optimize the embedding.

alldata <- FindNeighbors(alldata, reduction = "pca", graph.name = "SNN", assay = "RNA",
    k.param = 20, features = alldata@assays$RNA@var.features)

# Run UMAP on a graph

library(reticulate)
py_config()
#py_install("umap-learn")
#py_install("Numpy")
#use_python(python = "~/Anaconda3/bin/python", required = TRUE)
alldata <- RunUMAP(alldata, reduction.name = "UMAP_on_Graph", graph = "SNN", assay = "RNA")


#We can now plot the UMAP comparing both on PCA vs ScaledSata vs Graph.

p1 <- DimPlot(alldata, reduction = "umap", group.by = "orig.ident") + ggplot2::ggtitle(label = "UMAP_on_PCA")
p2 <- DimPlot(alldata, reduction = "UMAP_on_ScaleData", group.by = "orig.ident") +
    ggplot2::ggtitle(label = "UMAP_on_ScaleData")
p3 <- DimPlot(alldata, reduction = "UMAP_on_Graph", group.by = "orig.ident") + ggplot2::ggtitle(label = "UMAP_on_Graph")
leg <- get_legend(p1)

gridExtra::grid.arrange(gridExtra::arrangeGrob(p1 + NoLegend() + NoAxes(), p2 + NoLegend() +
    NoAxes(), p3 + NoLegend() + NoAxes(), leg, nrow = 2), ncol = 1, widths = c(1))




###Ploting genes of interest

myfeatures <- c("CD3E", "CD4", "CD8A", "NKG7", "GNLY", "MS4A1", "CD14", "LYZ", "MS4A7",
    "FCGR3A", "CST3", "FCER1A")


plot_list <- list()
for (i in myfeatures) {
    plot_list[[i]] <- FeaturePlot(alldata, reduction = "umap", dims = 1:2, features = i,
        ncol = 3, order = T) + NoLegend() + NoAxes() + NoGrid()
}
plot_grid(ncol = 3, plotlist = plot_list)

saveRDS(alldata, "data/results/seurat_qc_dr.rds")


####
###
#######
###integrating multiple single cell RNA-seq datasets


suppressPackageStartupMessages({
    library(Seurat)
    library(cowplot)
    library(ggplot2)
})

setwd("E:\\RNAseq\\scRNA\\pbmc_data\\pbmc_merge")

alldata <- readRDS("data/results/seurat_qc_dr.rds")
print(names(alldata@reductions))


#We split the combined object into a list, with each dataset as an element. 
#We perform standard preprocessing (log-normalization), 
#and identify variable features individually for 
#each dataset based on a variance stabilizing transformation (¡°vst¡±).


alldata.list <- SplitObject(alldata, split.by = "orig.ident")

for (i in 1:length(alldata.list)) {
    alldata.list[[i]] <- NormalizeData(alldata.list[[i]], verbose = FALSE)
    alldata.list[[i]] <- FindVariableFeatures(alldata.list[[i]], selection.method = "vst",
        nfeatures = 2000, verbose = FALSE)
}

hvgs_per_dataset <- lapply(alldata.list, function(x) {
    x@assays$RNA@var.features
})
# venn::venn(hvgs_per_dataset,opacity = .4,zcolor = scales::hue_pal()(3),cexsn
# = 1,cexil = 1,lwd=1,col='white',frame=F,borders = NA)

temp <- unique(unlist(hvgs_per_dataset))
overlap <- sapply(hvgs_per_dataset, function(x) {
    temp %in% x
})
pheatmap::pheatmap(t(overlap * 1), cluster_rows = F, color = c("grey90", "grey20"))

#
#We identify anchors using the FindIntegrationAnchors function, 
#which takes a list of Seurat objects as input.

alldata.anchors <- FindIntegrationAnchors(object.list = alldata.list, dims = 1:30,anchor.features = 5000,
    reduction = "cca")


#We then pass these anchors to the IntegrateData function, which returns a Seurat object.

alldata.int <- IntegrateData(anchorset = alldata.anchors, dims = 1:30, new.assay.name = "CCA")


#We can observe that a new assay slot is now created under the name CCA.


names(alldata.int@assays)

# by default, Seurat now sets the integrated assay as the default assay, so any
# operation you now perform will be on the ingegrated data.

alldata.int@active.assay


#After running IntegrateData, the Seurat object will contain a new Assay with the integrated 
#(or ¡®batch-corrected¡¯) expression matrix. 
#Note that the original (uncorrected values) are still stored in the object in the ¡°RNA¡± assay, 
#so you can switch back and forth. 
#We can then use this new integrated matrix for downstream analysis and visualization. 
#Here we scale the integrated data, run PCA, and visualize the results with UMAP and TSNE. 
#The integrated datasets cluster by cell type, instead of by technology.

# Run Dimensionality reduction on integrated space
alldata.int <- ScaleData(alldata.int, verbose = FALSE)
alldata.int <- RunPCA(alldata.int, npcs = 30, verbose = FALSE)
alldata.int <- RunUMAP(alldata.int, dims = 1:30)


alldata.int <- RunTSNE(alldata.int, dims = 1:30)


#We can now plot the un-integrated and the integrated space reduced dimensions.

plot_grid(ncol = 3,
  DimPlot(alldata, reduction = "pca", group.by = "orig.ident")+NoAxes()+ggtitle("PCA raw_data"),
  DimPlot(alldata, reduction = "tsne", group.by = "orig.ident")+NoAxes()+ggtitle("tSNE raw_data"),
  DimPlot(alldata, reduction = "umap", group.by = "orig.ident")+NoAxes()+ggtitle("UMAP raw_data"),
  
  DimPlot(alldata.int, reduction = "pca", group.by = "orig.ident")+NoAxes()+ggtitle("PCA integrated"),
  DimPlot(alldata.int, reduction = "tsne", group.by = "orig.ident")+NoAxes()+ggtitle("tSNE integrated"),
  DimPlot(alldata.int, reduction = "umap", group.by = "orig.ident")+NoAxes()+ggtitle("UMAP integrated")
)


plot_grid(ncol = 2,
#  DimPlot(alldata, reduction = "pca", group.by = "orig.ident")+NoAxes()+ggtitle("PCA raw_data"),
  DimPlot(alldata, reduction = "tsne", group.by = "orig.ident")+NoAxes()+ggtitle("tSNE raw_data"),
  DimPlot(alldata, reduction = "umap", group.by = "orig.ident")+NoAxes()+ggtitle("UMAP raw_data"),
  
#  DimPlot(alldata.int, reduction = "pca", group.by = "orig.ident")+NoAxes()+ggtitle("PCA integrated"),
  DimPlot(alldata.int, reduction = "tsne", group.by = "orig.ident")+NoAxes()+ggtitle("tSNE integrated"),
  DimPlot(alldata.int, reduction = "umap", group.by = "orig.ident")+NoAxes()+ggtitle("UMAP integrated")
)

###
###
###
#Let¡¯s plot some marker genes for different celltypes onto the embedding. Some genes are:

#Markers	Cell Type
#CD3E	T cells
#CD3E CD4	CD4+ T cells
#CD3E CD8A	CD8+ T cells
#GNLY, NKG7	NK cells
#MS4A1	B cells
#CD14, LYZ, CST3, MS4A7	CD14+ Monocytes
#FCGR3A, LYZ, CST3, MS4A7	FCGR3A+ Monocytes
#FCER1A, CST3	DCs

myfeatures <- c("CD3E", "CD4", "CD8A", "NKG7", "GNLY", "MS4A1", "CD14", "LYZ", "MS4A7",
    "FCGR3A", "CST3", "FCER1A")
plot_list <- list()
for (i in myfeatures) {
    plot_list[[i]] <- FeaturePlot(alldata.int, reduction = "umap", dims = 1:2, features = i,
        ncol = 3, order = T) + NoLegend() + NoAxes() + NoGrid()
}
plot_grid(ncol = 3, plotlist = plot_list)


###
####
####harmony

library(harmony)

alldata.harmony <- RunHarmony(alldata, group.by.vars = "orig.ident", reduction = "pca",
    dims.use = 1:50, assay.use = "RNA")

# Here we use all PCs computed from Harmony for UMAP calculation
alldata.int[["harmony"]] <- alldata.harmony[["harmony"]]
alldata.int <- RunUMAP(alldata.int, dims = 1:50, reduction = "harmony", reduction.name = "umap_harmony")
alldata.int <- RunTSNE(alldata.int, reduction = "harmony", reduction.name = "tsne_harmony")


hvgs <- unique(unlist(hvgs_per_dataset))

assaylist <- list()
genelist <- list()
for (i in 1:length(alldata.list)) {
    assaylist[[i]] <- t(as.matrix(GetAssayData(alldata.list[[i]], "data")[hvgs, ]))
    genelist[[i]] <- hvgs
}

lapply(assaylist, dim)


####
####
####scanorama

library(reticulate)
py_config()
#pip install scanorama
#reticulate::py_install('scanorama', pip = T)

scanorama <- import("scanorama")

integrated.data <- scanorama$integrate(datasets_full = assaylist, genes_list = genelist)

intdimred <- do.call(rbind, integrated.data[[1]])
colnames(intdimred) <- paste0("PC_", 1:100)
rownames(intdimred) <- colnames(alldata.int)

# Add standard deviations in order to draw Elbow Plots in Seurat
stdevs <- apply(intdimred, MARGIN = 2, FUN = sd)

alldata.int[["scanorama"]] <- CreateDimReducObject(embeddings = intdimred, stdev = stdevs,
    key = "PC_", assay = "RNA")

#Here we use all PCs computed from Scanorama for UMAP calculation
alldata.int <- RunUMAP(alldata.int, dims = 1:100, reduction = "scanorama", reduction.name = "umap_scanorama")
alldata.int <- RunTSNE(alldata.int, reduction = "scanorama", reduction.name = "tsne_scanorama")


p1 <- DimPlot(alldata, reduction = "umap", group.by = "orig.ident") + ggtitle("UMAP raw_data")
p2 <- DimPlot(alldata.int, reduction = "umap", group.by = "orig.ident") + ggtitle("UMAP CCA")
p3 <- DimPlot(alldata.int, reduction = "umap_harmony", group.by = "orig.ident") +
    ggtitle("UMAP Harmony")
p4 <- DimPlot(alldata.int, reduction = "umap_scanorama", group.by = "orig.ident") +
    ggtitle("UMAP Scanorama")
leg <- get_legend(p1)

gridExtra::grid.arrange(gridExtra::arrangeGrob(p1 + NoLegend() + NoAxes(), p2 + NoLegend() +
    NoAxes(), p3 + NoLegend() + NoAxes(), p4 + NoLegend() + NoAxes(), nrow = 2),
    leg, ncol = 2, widths = c(8, 2))



p1 <- DimPlot(alldata, reduction = "tsne", group.by = "orig.ident") + ggtitle("tsne raw_data")
p2 <- DimPlot(alldata.int, reduction = "tsne", group.by = "orig.ident") + ggtitle("tsne CCA")
p3 <- DimPlot(alldata.int, reduction = "tsne_harmony", group.by = "orig.ident") +
    ggtitle("tsne Harmony")
p4 <- DimPlot(alldata.int, reduction = "tsne_scanorama", group.by = "orig.ident") + ggtitle("tsne Scanorama")
leg <- get_legend(p1)

gridExtra::grid.arrange(gridExtra::arrangeGrob(p1 + NoLegend() + NoAxes(), 
    p2 + NoLegend() + NoAxes(), p3 + NoLegend() + NoAxes(), 
    p4 + NoLegend() + NoAxes(), 
    nrow = 2),
    leg, ncol = 2, widths = c(8, 2))




saveRDS(alldata.int, "data/results/seurat_qc_dr_int.rds")

#alldata.int <- readRDS("data/results/seurat_qc_dr_int.rds")
alldata.int
head(alldata.int@meta.data)
DefaultAssay(alldata.int)
as.data.frame(alldata.int@assays$RNA@counts[1:10, 1:5])
as.data.frame(alldata.int@assays$CCA[1:10, 1:5])
#as.data.frame(alldata.int@assays$CCA[1:10, 1:5])
print(names(alldata.int@reductions))
#as.data.frame(alldata.int@reductions$umap)
str(alldata.int)
DefaultAssay(alldata.int)


alldata
head(alldata@meta.data)
as.data.frame(alldata@assays$RNA@counts[1:10, 1:2])

print(names(alldata@reductions))
str(alldata)

DefaultAssay(alldata)





#################
############
###########Clustering

#Let¡¯s first load all necessary libraries and also the integrated dataset from the previous step.

if (!require(clustree)) {
    install.packages("clustree", dependencies = FALSE)
}

suppressPackageStartupMessages({
    library(Seurat)
    library(cowplot)
    library(ggplot2)
    library(pheatmap)
    library(rafalib)
    library(clustree)
})

setwd("E:\\RNAseq\\scRNA\\pbmc_data\\pbmc_merge")

alldata <- readRDS("data/results/seurat_qc_dr_int.rds")

alldata@assays$CCA
#The procedure of clustering on a Graph can be generalized as 3 main steps:
#Build a kNN graph from the data
#Prune spurious connections from kNN graph (optional step). This is a SNN graph.
#Find groups of cells that maximizes the connections within the group compared other groups.


#The first step into graph clustering is to construct a k-nn graph, 
#in case you don¡¯t have one. For this, we will use the PCA space. 
#Thus, as done for dimensionality reduction, 
#we will use ony the top N PCA dimensions for this purpose (the same used for computing UMAP / tSNE).

# check that CCA is still the active assay
alldata@active.assay

alldata <- FindNeighbors(alldata, dims = 1:30, k.param = 60, prune.SNN = 1/15)

names(alldata@graphs)


#We can take a look at the kNN graph. 
#It is a matrix where every connection between cells is represented as 1s. 

pheatmap(alldata@graphs$CCA_nn[1:200, 1:200], col = c("white", "black"), border_color = "grey90",
    legend = F, cluster_rows = F, cluster_cols = F, fontsize = 2)



#Clustering on a graph

#Once the graph is built, we can now perform graph clustering. 
#The clustering is done respective to a resolution which can be interpreted as how coarse you want your cluster to be. 
#Higher resolution means higher number of clusters.

#In Seurat, the function FindClusters will do a graph-based clustering using ¡°Louvain¡± algorithim by default (algorithm = 1). 
#TO use the leiden algorithm, you need to set it to algorithm = 4.


for (res in c(0.1,0.15,0.2, 0.25, 0.3,0.4,0.5,0.75, 1)) {
    alldata <- FindClusters(alldata, graph.name = "CCA_snn", resolution = res, algorithm = 1)
}

# each time you run clustering, the data is stored in meta data columns:
# seurat_clusters - lastest results only CCA_snn_res.XX - for each different
# resolution you test.


plot_grid(ncol = 3, DimPlot(alldata, reduction = "umap", group.by = "CCA_snn_res.0.1") + ggtitle("louvain_0.1"), 
                    DimPlot(alldata, reduction = "umap", group.by = "CCA_snn_res.0.15") + ggtitle("louvain_0.15"),
                    DimPlot(alldata, reduction = "umap", group.by = "CCA_snn_res.0.2") + ggtitle("louvain_0.2"), 
                    DimPlot(alldata, reduction = "umap", group.by = "CCA_snn_res.0.25") + ggtitle("louvain_0.25"),
                    DimPlot(alldata, reduction = "umap", group.by = "CCA_snn_res.0.3") + ggtitle("louvain_0.3"),
                    DimPlot(alldata, reduction = "umap", group.by = "CCA_snn_res.0.4") + ggtitle("louvain_0.4"),
DimPlot(alldata, reduction = "umap", group.by = "CCA_snn_res.0.5") + ggtitle("louvain_0.5"),
DimPlot(alldata, reduction = "umap", group.by = "CCA_snn_res.0.75") + ggtitle("louvain_0.75"),
DimPlot(alldata, reduction = "umap", group.by = "CCA_snn_res.1") + ggtitle("louvain_1")
         )


#We can now use the clustree package to visualize 
#how cells are distributed between clusters depending on resolution.

suppressPackageStartupMessages(library(clustree))

clustree(alldata@meta.data, prefix = "CCA_snn_res.")


###K-means clustering
###
###

for (k in c(5, 7, 10, 12, 15, 17, 20)) {
    alldata@meta.data[, paste0("kmeans_", k)] <- kmeans(x = alldata@reductions[["pca"]]@cell.embeddings,
        centers = k, nstart = 100)$cluster
}


plot_grid(ncol = 2, DimPlot(alldata, reduction = "umap", group.by = "kmeans_5") +
    ggtitle("kmeans_5"), DimPlot(alldata, reduction = "umap", group.by = "kmeans_7") +
    ggtitle("kmeans_7"), DimPlot(alldata, reduction = "umap", group.by = "kmeans_10") +
    ggtitle("kmeans_10"), DimPlot(alldata, reduction = "umap", group.by = "kmeans_15") +
    ggtitle("kmeans_15"))

clustree(alldata@meta.data, prefix = "kmeans_")



###Hierarchical clustering
####Defining distance between cells

d <- dist(alldata@reductions[["pca"]]@cell.embeddings, method = "euclidean")

# Compute sample correlations
sample_cor <- cor(Matrix::t(alldata@reductions[["pca"]]@cell.embeddings))

# Transform the scale from correlations
sample_cor <- (1 - sample_cor)/2

# Convert it to a distance object
d2 <- as.dist(sample_cor)

# euclidean
h_euclidean <- hclust(d, method = "ward.D2")

# correlation
h_correlation <- hclust(d2, method = "ward.D2")

#euclidean distance
alldata$hc_euclidean_5 <- cutree(h_euclidean,k = 5)
alldata$hc_euclidean_10 <- cutree(h_euclidean,k = 10)
alldata$hc_euclidean_15 <- cutree(h_euclidean,k = 15)

#correlation distance
alldata$hc_corelation_5 <- cutree(h_correlation,k = 5)
alldata$hc_corelation_10 <- cutree(h_correlation,k = 10)
alldata$hc_corelation_15 <- cutree(h_correlation,k = 15)


plot_grid(ncol = 3,
  DimPlot(alldata, reduction = "umap", group.by = "hc_euclidean_5")+ggtitle("hc_euc_5"),
  DimPlot(alldata, reduction = "umap", group.by = "hc_euclidean_10")+ggtitle("hc_euc_10"),
  DimPlot(alldata, reduction = "umap", group.by = "hc_euclidean_15")+ggtitle("hc_euc_15"),

  DimPlot(alldata, reduction = "umap", group.by = "hc_corelation_5")+ggtitle("hc_cor_5"),
  DimPlot(alldata, reduction = "umap", group.by = "hc_corelation_10")+ggtitle("hc_cor_10"),
  DimPlot(alldata, reduction = "umap", group.by = "hc_corelation_15")+ggtitle("hc_cor_15")
)


head(alldata@meta.data)

saveRDS(alldata, "data/results/seurat_qc_dr_int_cl.rds")



###
#Differential gene expression
####
####
####

suppressPackageStartupMessages({
  library(Seurat)
  library(venn)
  library(dplyr)
  library(cowplot)
  library(ggplot2)
  library(pheatmap)
  library(enrichR)
  library(rafalib)
})

setwd("E:\\RNAseq\\scRNA\\pbmc_data\\pbmc_merge")

alldata <- readRDS("data/results/seurat_qc_dr_int_cl.rds")

####select by annotaion
head(alldata@meta.data)
names(alldata@meta.data)

print(names(alldata@reductions))

sel.clust = "CCA_snn_res.0.4"

plot_grid(ncol = 2,
  DimPlot(alldata, reduction = "tsne", group.by = "orig.ident") +ggtitle ("orig.ident") ,
  DimPlot(alldata, reduction = "tsne", label = T, group.by = sel.clust)  + ggtitle(sel.clust),
  DimPlot(alldata, reduction = "tsne_harmony", label = T, group.by = sel.clust)  + ggtitle(sel.clust),
   DimPlot(alldata, reduction = "tsne_scanorama", label = T, group.by = sel.clust)  + ggtitle(sel.clust)
  )


#meta <- alldata@meta.data

#clu_list <- split(meta$cell_annotation, meta[,sel.clust])
#lapply(clu_list,length) 
#lapply(clu_list,function(x) {
#                 table(x) %>% sort(decreasing = T)}) 


#Set the identity as louvain with resolution 0.25
#sel.clust = "CCA_snn_res.0.4"

alldata <- SetIdent(alldata, value = sel.clust)
table(alldata@active.ident)

# plot this clustering
plot_grid(ncol = 1,
  DimPlot(alldata, label = T) + NoAxes(),
  DimPlot(alldata, group.by = "orig.ident") + NoAxes()
  ) 

#Cell marker genes
#Let us first compute a ranking for the highly differential genes in each cluster. 
#There are many different tests and parameters to be chosen that can be used to refine your results. 
#When looking for marker genes, we want genes that are positivelly expressed in 
#a cell type and possibly not expressed in the others.

#Compute differentiall expression
#markers <- FindAllMarkers(object, test.use = 'LR', latent.vars = 'batch') ##I think batch use LR is better




markers_genes <- FindAllMarkers(alldata,
#                               logfc.threshold = 0.2,
                               test.use = "wilcox",
                               min.pct = 0.1,
#                               min.diff.pct = 0.2,
#                               only.pos = TRUE,
                               max.cells.per.ident = 50,
#                               assay = "RNA")
                               assay = "CCA")



library(openxlsx)
gwb <- createWorkbook()
addWorksheet(gwb,sheetName = 'single_cell_mark')
writeData(gwb,sheet = 'single_cell_mark',markers_genes)

saveWorkbook(gwb, "DEGs.xlsx", overwrite = TRUE)


#markers_genes <- FindAllMarkers(alldata,
#                              logfc.threshold = 0.2,
#                               test.use = "LR",
#                               latent.vars = 'orig.ident',
#                               min.pct = 0.1,
#                               min.diff.pct = 0.2,
#                               only.pos = TRUE,
#                               max.cells.per.ident = 50,
 #                              assay = "RNA")

markers_genes %>% filter(avg_log2FC > 0) %>%
                  group_by(cluster)  %>% top_n(-25, p_val_adj) -> top25
top25

#We can now select the top 25 up regulated genes for plotting.
mypar(2,5,mar=c(4,6,3,1))
for(i in unique(top25$cluster)){
  barplot( sort( setNames(top25$avg_log2FC, top25$gene) [top25$cluster == i], F),
           horiz = T,las=1 ,main=paste0(i," vs. rest"),border = "white", yaxs="i" )
  abline(v=c(0,0.25),lty=c(1,2))
}

####We can visualize them as a heatmap. Here we are selecting the top 5.
markers_genes %>% filter(avg_log2FC > 0) %>% group_by(cluster)  %>% 
                  top_n(-5, p_val_adj) -> top5

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

saveRDS(alldata,"data/results/seurat_qc_dr_int_cl_dge")

##
##Differential expression across conditions

#The second way of computing differential expression is to answer which genes are differentially expressed within a cluster. For example, in our case we have libraries comming from patients and controls and we would like to know which genes are influenced the most in a particular cell type.

#For this end, we will first subset our data for the desired cell cluster, 
#then change the cell identities to the variable of comparison 
#(which now in our case is the ¡°type¡±, e.g. Covid/Ctrl).
# select all cells in cluster 1
#cell_selection <- subset(alldata, cells = colnames(alldata)[ alldata@meta.data[,sel.clust] == 2])
#cell_selection <- SetIdent(cell_selection, value = "type")
#Compute differentiall expression
#DGE_cell_selection <- FindAllMarkers(cell_selection,
#                               logfc.threshold = 0.2,
#                               test.use = "wilcox",
#                               min.pct = 0.1,
#                               min.diff.pct = 0.2,
#                               only.pos = TRUE,
#                               max.cells.per.ident = 50,
#                               assay = "RNA")


#We can now plot the expression across the ¡°type¡±.
#DGE_cell_selection %>% group_by(cluster)  %>% top_n(-5, p_val) -> top5_cell_selection

#VlnPlot(cell_selection, features = as.character(unique(top5_cell_selection$gene)),
#        ncol = 5,group.by = "type",assay = "RNA", pt.size = .1)

#We can also plot these genes across all clusters, 
#but split by ¡°type¡±, to check if the genes are also up/downregulated in other celltypes.

#VlnPlot(alldata, features = as.character(unique(top5_cell_selection$gene)),
#        ncol = 5, split.by = "type",assay = "RNA", pt.size = 0)


#Gene Set Analysis
# Load additional packages
library(enrichR)

# Check available databases to perform enrichment (then choose one)
db <- enrichR::listEnrichrDbs()
sort(db$libraryName)
# Perform enrichment
genes <- markers_genes %>% filter(avg_log2FC > 0) %>% filter(cluster == 0 & p_val_adj < 0.05) 
enrich_results <- enrichr(
 genes     =  genes$gene,
 databases = "GO_Biological_Process_2021" )[[1]]


#Gene Set Enrichment Analysis (GSEA)

#gene_rank <- setNames( DGE_cell_selection$avg_log2FC, casefold(rownames(DGE_cell_selection),upper=T) )

gene_rank <- setNames(genes$avg_log2FC, casefold(rownames(genes),upper=T) )

# install.packages("msigdbr")
library(msigdbr)

#Download gene sets
msigdbgmt <- msigdbr::msigdbr("Homo sapiens")
msigdbgmt <- as.data.frame(msigdbgmt)

#List available gene sets
unique(msigdbgmt$gs_subcat)

#Subset which gene set you want to use.
msigdbgmt_subset <- msigdbgmt[msigdbgmt$gs_subcat == "CP:KEGG",]
gmt <- lapply( unique(msigdbgmt_subset$gs_name),function(x){msigdbgmt_subset [msigdbgmt_subset$gs_name == x ,"gene_symbol"]} )
names(gmt) <- unique(paste0(msigdbgmt_subset$gs_name,"_",msigdbgmt_subset$gs_exact_source))

#Next, we will be using the GSEA. This will result in a table containing information for several pathways. 
#We can then sort and filter those pathways to visualize only the top ones. 
#You can select/filter them by either p-value or normalized enrichemnet score (NES).

#saveRDS(alldata,"data/results/seurat_qc_dr_int_cl_dge")
#write.csv(markers_genes)

#saveRDS(alldata,"data/results/seurat_qc_dr_int_cl_dge")
###
##########cell type gsea


suppressPackageStartupMessages({
    library(dplyr)
    library(cowplot)
    library(rafalib)
    library(scPred)
})


alldata
head(alldata@meta.data)
table(alldata@meta.data["orig.ident"])
table(alldata$orig.ident)
as.data.frame(alldata@assays$RNA@counts[1:10, 1:2])
print(names(alldata@reductions))
#str(alldata)
DefaultAssay(alldata)

slice <- 1:30

#DGE_list <- markers_genes %>% filter(p_val < 0.001) %>% group_by(cluster) %>%top_n (-80,p_val)

DGE_list <- markers_genes %>% filter(!is.na(gene)) %>%
            filter(avg_log2FC > 0) %>% arrange(p_val) %>%
            filter(p_val < 0.001) %>% group_by(cluster) %>% slice(slice)


DGE_list <- split(DGE_list, DGE_list$cluster)
unlist(lapply(DGE_list, nrow))


# Load the human marker table
library(openxlsx)
markers <- read.xlsx("E:\\RNAseq\\scRNA\\db\\CellMarker_2\\Cell_marker_Human.xlsx",sheet=1)

markers <- markers[markers$species == "Human", ]
markers <- markers[markers$cell_type == "Normal cell", ]


# Filter by tissue (to reduce computational time and have tissue-specific
# classification) sort(unique(markers$tissueType))
# grep('blood',unique(markers$tissueType),value = T) markers <- markers [
# markers$tissueType %in% c('Blood','Venous blood', 'Serum','Plasma',
# 'Spleen','Bone marrow','Lymph node'), ]

# remove strange characters etc.
celltype_list <- lapply(unique(markers$cell_name), function(x) {
    x <- paste(markers$Symbol[markers$cell_name == x], sep = ",")
    x <- gsub("[[]|[]]| |-", ",", x)
    x <- unlist(strsplit(x, split = ","))
    x <- unique(x[!x %in% c("", "NA", "family")])
    x <- casefold(x, upper = T)
})

names(celltype_list) <- unique(markers$cell_name)
# celltype_list <- lapply(celltype_list , function(x) {x[1:min(length(x),50)]}
# )

#celltype_list <- celltype_list[unlist(lapply(celltype_list, length)) < 100]
#celltype_list <- celltype_list[unlist(lapply(celltype_list, length)) > 3]



library(fgsea)
# run fgsea for each of the clusters in the list
res <- lapply(DGE_list, function(x) {
    gene_rank <- setNames(x$avg_log2FC, x$gene)
#   maxSize = length(stats) - 1,    
    fgseaRes <- fgsea(pathways = celltype_list, stats = gene_rank,
    maxSize = length(gene_rank) - 1, nperm = 10000)
    return(fgseaRes)
})
names(res) <- names(DGE_list)


# You can filter and resort the table based on ES, NES or pvalue
#res <- lapply(res, function(x) {
#    x[x$pval < 0.1, ]
#})
res <- lapply(res, function(x) {
    x[x$NES > 0, ]
})


res <- lapply(res, function(x) {
    x[order(x$pval), ]
})

res <- lapply(res, function(x) {
    x[order(x$size, decreasing = T), ]
})

res <- lapply(res, function(x) {
    x[x$size > 1, ]
})

# show top 3 for each cluster.
lapply(res, head, 3)

cluster_n <- lapply(res,nrow) %>% unlist

cell_marker_gsea <- dplyr::bind_rows(res, .id = "cluster")  ###change dataframe , add column cluster by list name

library(openxlsx)
gwb <- createWorkbook("")

markers <- markers_genes %>% filter(avg_log2FC > 0) %>% group_by(cluster) %>% slice(slice)


addWorksheet(gwb,sheetName = 'cell_marker')
writeData(gwb,sheet = 'cell_marker',markers)

addWorksheet(gwb,sheetName = 'cell_marker_gsea')
writeData(gwb,sheet = 'cell_marker_gsea',cell_marker_gsea)

saveWorkbook(gwb, "cell_marker_gsea.xlsx", overwrite = TRUE)


### pbmc paper

slice <- 1:100

library(openxlsx)
gwb <- createWorkbook("")

markers <- markers_genes %>% filter(avg_log2FC > 0) %>% group_by(cluster) %>% slice(slice)

addWorksheet(gwb,sheetName = 'cell_marker')
writeData(gwb,sheet = 'cell_marker',markers)

#DGE_list <- markers_genes %>% filter(p_val < 0.001) %>% group_by(cluster) %>%top_n (-80,p_val)

DGE_list <- markers_genes %>% filter(!is.na(gene)) %>%
            filter(avg_log2FC > 0) %>% arrange(p_val) %>%
            filter(p_val < 0.001) %>% group_by(cluster) %>% slice(slice)


DGE_list <- split(DGE_list, DGE_list$cluster)
unlist(lapply(DGE_list, nrow))


markers <- read.xlsx("NIHMS1585672-supplement-Supp_Supp__Tables_1-12_and_15-19.xlsx",sheet="Table s10")
markers_p <- markers %>% filter(Positive.or.Negative == "+")
markers_n <- markers %>% filter(Positive.or.Negative == "-")

# remove strange characters etc.
celltype_list_p <- lapply(unique(markers_p$Cell.type), function(x) {
    x <- paste(markers_p$Gene[markers_p$Cell.type == x], sep = ",")
#    x <- gsub("[[]|[]]| |-", ",", x)
    x <- unlist(strsplit(x, split = ","))
#    x <- unique(x[!x %in% c("", "NA", "family")])
    x <- casefold(x, upper = T)
})

names(celltype_list_p) <- unique(markers_p$Cell.type)

####


library(fgsea)
# run fgsea for each of the clusters in the list
res <- lapply(DGE_list, function(x) {
    gene_rank <- setNames(x$avg_log2FC, x$gene)
#   maxSize = length(stats) - 1,    
    fgseaRes <- fgsea(pathways = celltype_list_p, stats = gene_rank,
    maxSize = length(gene_rank) - 1, nperm = 10000)
    return(fgseaRes)
})
names(res) <- names(DGE_list)


# You can filter and resort the table based on ES, NES or pvalue
res <- lapply(res, function(x) {
    x[x$pval < 0.1, ]
})
res <- lapply(res, function(x) {
    x[x$NES > 0, ]
})


res <- lapply(res, function(x) {
    x[order(x$pval), ]
})

res <- lapply(res, function(x) {
    x[order(x$size, decreasing = T), ]
})

res <- lapply(res, function(x) {
    x[x$size > 1, ]
})

# show top 3 for each cluster.
lapply(res, head, 3)

cluster_n <- lapply(res,nrow) %>% unlist

cell_marker_gsea <- dplyr::bind_rows(res, .id = "cluster")  ###change dataframe , add column cluster by list name



addWorksheet(gwb,sheetName = 'cell_marker_p')
writeData(gwb,sheet = 'cell_marker_p',cell_marker_gsea)



# remove strange characters etc.
celltype_list_n <- lapply(unique(markers_n$Cell.type), function(x) {
    x <- paste(markers_n$Gene[markers_n$Cell.type == x], sep = ",")
#    x <- gsub("[[]|[]]| |-", ",", x)
    x <- unlist(strsplit(x, split = ","))
#    x <- unique(x[!x %in% c("", "NA", "family")])
    x <- casefold(x, upper = T)
})

names(celltype_list_n) <- unique(markers_n$Cell.type)

####


library(fgsea)
# run fgsea for each of the clusters in the list
res <- lapply(DGE_list, function(x) {
    gene_rank <- setNames(x$avg_log2FC, x$gene)
#   maxSize = length(stats) - 1,    
    fgseaRes <- fgsea(pathways = celltype_list_n, stats = gene_rank,
    maxSize = length(gene_rank) - 1, nperm = 10000)
    return(fgseaRes)
})
names(res) <- names(DGE_list)


# You can filter and resort the table based on ES, NES or pvalue
res <- lapply(res, function(x) {
    x[x$pval < 0.1, ]
})
res <- lapply(res, function(x) {
    x[x$NES > 0, ]
})


res <- lapply(res, function(x) {
    x[order(x$pval), ]
})

res <- lapply(res, function(x) {
    x[order(x$size, decreasing = T), ]
})

res <- lapply(res, function(x) {
    x[x$size > 1, ]
})

# show top 3 for each cluster.
lapply(res, head, 3)

cluster_n <- lapply(res,nrow) %>% unlist

cell_marker_gsea <- dplyr::bind_rows(res, .id = "cluster")  ###change dataframe , add column cluster by list name



addWorksheet(gwb,sheetName = 'cell_marker_n')
writeData(gwb,sheet = 'cell_marker_n',cell_marker_gsea)


##meatadata

meta_cell_type <- read.xlsx("single_cell_metadata_2000.xlsx",sheet=1) %>% select(cell_id, cell_type_final)

meta_df <- alldata@meta.data %>% mutate(cell_id = rownames(.)) %>% left_join(.,meta_cell_type) %>% select(sel.clust, cell_type_final)

meta_df <- split(meta_df[,2], meta_df[,1]) %>%

                  lapply (function(x) {
                                     x <- data.frame(table(x))
                                     names(x)[1] <- "cell_type_2000"
                                     x <- arrange(x,desc(Freq))
                  }
                  ) %>%
                  dplyr::bind_rows(., .id = sel.clust)

addWorksheet(gwb,sheetName = 'cell_marker_2000')
writeData(gwb,sheet = 'cell_marker_2000',meta_df)

saveWorkbook(gwb, "cell_marker_paper.xlsx", overwrite = TRUE)


plot_grid(ncol = 2,
#  DimPlot(alldata, reduction = "tsne", group.by = "orig.ident") +ggtitle ("orig.ident") ,
  DimPlot(alldata, reduction = "tsne", label = T, group.by = sel.clust)  + ggtitle(sel.clust),
  DimPlot(alldata, reduction = "tsne_harmony", label = T, group.by = sel.clust)  + ggtitle(sel.clust),
   DimPlot(alldata, reduction = "tsne_scanorama", label = T, group.by = sel.clust)  + ggtitle(sel.clust)
  )


DimPlot(alldata, reduction = "tsne", label = T, group.by = sel.clust) +ggtitle ("orig.ident")

unique(markers$Cell.type)




##########cell type annotantion

############do literture mining by gsea result and then come back
######
#####
##########cell type prediction by mining results
########
#Celltype prediction

#Celltype prediction can either be performed on 
#indiviudal cells where each cell gets a predicted celltype label, 
#or on the level of clusters. All methods are based on similarity to other datasets, 
#single cell or sorted bulk RNAseq, or uses know marker genes for each celltype.


suppressPackageStartupMessages({
    library(Seurat)
    library(dplyr)
    library(cowplot)
    library(ggplot2)
    library(pheatmap)
    library(rafalib)
    library(scPred)
    library(openxlsx)
})

setwd("E:\\RNAseq\\scRNA\\pbmc_data\\pbmc_merge")

alldata <- readRDS("data/results/seurat_qc_dr_int_cl_dge")
### data information
alldata
head(alldata@meta.data)
table(alldata@meta.data["orig.ident"])
table(alldata$orig.ident)
as.data.frame(alldata@assays$RNA@counts[1:10, 1:2])
print(names(alldata@reductions))
#str(alldata)
DefaultAssay(alldata)
###
sel.clust = "CCA_snn_res.0.4"
alldata$sel.clust <- alldata$CCA_snn_res.0.4



plot_grid(ncol = 1,
  DimPlot(alldata, reduction = "tsne", group.by = "orig.ident") +ggtitle ("Site") ,
  DimPlot(alldata, reduction = "tsne", label = T, group.by = sel.clust)  + ggtitle(sel.clust)
  )


#meta <- alldata@meta.data
#clu_list <- split(meta$cell_annotation, meta[,sel.clust])
#lapply(clu_list,length) 
#lapply(clu_list,function(x) {
#                 table(x) %>% sort(decreasing = T)}) 

#obj <- createscCATCH(data = Seurat_obj[['RNA']]@data, cluster = as.character(Idents(Seurat_obj)))




cowplot::plot_grid(ncol = 2, 
#                             DimPlot(alldata, label = T, reduction= "tsne",group.by = "sel.clust") + NoAxes() + NoLegend() + ggtitle("cluster_tnse"),
                             DimPlot(alldata, label = T, reduction= "umap",group.by = "sel.clust") + NoAxes() + NoLegend() + ggtitle("cluster_umap"),
                             DimPlot(alldata, label = T, reduction= "umap_harmony",group.by = "sel.clust") + NoAxes() + NoLegend() + ggtitle("cluster_harmony"),
                             DimPlot(alldata, label = T, reduction= "umap_scanorama",group.by = "sel.clust") + NoAxes() + NoLegend() + ggtitle("cluster_scanorama")
                 )

cowplot::plot_grid(ncol = 2, 
                             DimPlot(alldata, label = T, reduction= "tsne",group.by = "sel.clust") + NoAxes() + NoLegend() + ggtitle("cluster_tnse"),
 #                            DimPlot(alldata, label = T, reduction= "umap",group.by = "sel.clust") + NoAxes() + NoLegend() + ggtitle("cluster_umap"),
                             DimPlot(alldata, label = T, reduction= "tsne_harmony",group.by = "sel.clust") + NoAxes() + NoLegend() + ggtitle("cluster_harmony"),
                             DimPlot(alldata, label = T, reduction= "tsne_scanorama",group.by = "sel.clust") + NoAxes() + NoLegend() + ggtitle("cluster_scanorama")
                 )



fplot <- function(gene){
    FeaturePlot(alldata, reduction = "tsne", dims = 1:2, features = gene, min.cutoff = 0.1,max.cutoff = 8,
#         blend = TRUE,
#        ncol = 1, order = T) + NoLegend() + NoAxes() + NoGrid() +
         ncol = 1, order = T)  + NoAxes() + NoGrid() +
#        scale_color_gradientn( colours = c('lightgrey', 'red'), na.value = "lightgrey",  limits = limit)
#        scale_colour_gradientn(colours = rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")))
        scale_colour_gradientn(colours = colorRampPalette(colors = c("grey","orange","red"))(11))
}

plot_grid(ncol = 1,
#  DimPlot(alldata, reduction = "tsne", group.by = "orig.ident") +ggtitle ("Site") ,
  DimPlot(alldata, reduction = "tsne", label = T, group.by = sel.clust)  + ggtitle(sel.clust),
  fplot("IL2RA") ##CD25 
  )
fplot("TGFB1")


fplot("CD1D")

fplot("CD69")

fplot("CXCR3")

fplot("NOD1")


plot_list <- list()
plot_list[[1]] <- fplot("CD4")
plot_list[[2]] <- fplot("CD8A")
plot_list[[3]] <- fplot("CD8B")
plot_grid(ncol = 1, plotlist = plot_list)


fplot("CD4") 

fplot("CD3G") 
fplot("CD3D") 

fplot("CD3E") 


fplot("NCAM1")  ##CD56

fplot("FCGR3A") ##CD16

fplot("CXCR3")
fplot("CCR6")
fplot("KLRB1") ##CD161
fplot("IL2RA") ##CD25 

fplot("IL10") 
fplot("TGFB1") 


fplot("CXCR5")

fplot("FCGR3A") ##CD16
fplot("CD14")
fplot("CXCR4")


cell_type <- read.xlsx("cell_annotation.xlsx",sheet=1)


#alldata@meta.data <- left_join(alldata@meta.data,cell_type)

##########
#new.cluster.ids <- unlist(lapply(res, function(x) {
#    as.data.frame(x)[1, 1]
#})rownam)

dim(cell_type)
rownames(cell_type) <- cell_type[,1]

cell_type <- cell_type[alldata@meta.data$sel.clust,]

alldata@meta.data <- cbind(alldata@meta.data,cell_type)

#alldata@meta.data <- new.cluster.ids[as.character(alldata@active.ident)]

#alldata@reduction = "pca"

#names(alldata@reductions)


cowplot::plot_grid(ncol = 1, DimPlot(alldata, label = T, reduction= "umap",group.by = "pathway") + NoAxes() + NoLegend() + ggtitle("cluster"),
                             DimPlot(alldata, label = T, reduction= "tsne",group.by = "pathway") + NoAxes() + NoLegend() + ggtitle("cell_type"))

cowplot::plot_grid(ncol = 1, DimPlot(alldata, label = T, reduction= "tsne",group.by = "sel.clust") + NoAxes() + NoLegend() + ggtitle("cluster"),
                             DimPlot(alldata, label = T, reduction= "tsne",group.by = "pathway") + NoAxes() + NoLegend() + ggtitle("cell_type"))



#DimPlot(alldata, label = T, label.size = 4,reduction= "umap_harmony",group.by = "pathway") + NoAxes() + NoLegend() + ggtitle("cell_type") 
#labels=function(x) str_wrap(x, width=50)

#alldata@meta.data$pathway <- gsub("Plasmacytoid dendritic cell.*","pDC",alldata@meta.data$pathway)
#alldata@meta.data$pathway <- gsub("Natural killer T.*","NKT cell",alldata@meta.data$pathway)
table(alldata@meta.data$pathway )
#alldata <- SetIdent(alldata, value = sel.clust)

alldata@meta.data$cell_type_final <- alldata@meta.data$pathway

#DimPlot(alldata, label = T, label.size = 3.8,reduction= "umap_harmony",group.by = "cell_type_final") + NoAxes() + NoLegend() + ggtitle("cell_type") 

saveRDS(alldata,"data/results/seurat_qc_dr_int_cl_dge_cell_type")

gwb <- createWorkbook("")


addWorksheet(gwb,sheetName = 'metadata')
writeData(gwb,sheet = 'metadata',alldata@meta.data)



saveWorkbook(gwb, "single_cell_metadata.xlsx", overwrite = TRUE)






























