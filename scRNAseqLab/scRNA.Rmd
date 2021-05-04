---
title: "scRNA_Lab"
author: "Sid"
date: "5/3/2021"
output:
  html_document: default
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_chunk$set(fig.width=12, fig.height=8) 
```

# Introduction 
This week's lab was to follow the **Seurat Clustering Tutorial**. This tutorial uses Peripheral Blood Mononuclear Cell (PBMC) data provided by 10X Genomics. The following code is borrowed from the Seurat tutorial which can be found here: <https://satijalab.org/seurat/articles/pbmc3k_tutorial.html>. 

```{r, message=FALSE} 
## If libraries are not pre-installed then uncomment the lines below to install them 
# install.packages('dplyr')
# install.packages('Seurat')
# install.packages('patchwork')

# import necessary libraries 
library(dplyr)
library(Seurat)
library(patchwork)

# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "~/Documents/2020-21/BCB 5250/filtered_gene_bc_matrices/hg19/")
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc
```

# Preprocessing: Quality Control 
As stated by the authors, the following steps are the standard preprocessing measures taken in by Seurat. First, for quality control (QC) I filtered cells that had the following feature counts (FCs) and mitochondrial counts (mtCs): 
  
  1. FCs < 200
  2. FCs > 2,500
  3. mtCs > 5% 

```{r, message=FALSE} 
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```

Above we can see violin plots showing the frequencies and distributions of FCs, number of counts, and mtCs. 
```{r, message=FALSE} 
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
```

I applied the further QC which the authors recommended to visualize the inter-feature relationships of the data. The left-hand plot shows the relationship between mtCs and number of RNAs, while the latter shows the relationship between the FCs and number of RNAs. Below I applied the first filter to the data to remove observations based on the criteria provided above and from the visualizations. 
```{r, message=FALSE} 
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
```

# Preprocessing: Normalization 

Normalization is a common step of data preprocessing as it often speeds up computational run times. The authors suggested a global-scaling method, *log normalize* which used the total expression multiplied by a scalar (10,000 if default param) to normalize the feature expression for each cell before log-transforming the resultant value. I followed the author's code exactly. 

```{r, message=FALSE} 
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- NormalizeData(pbmc)
```

# Preprocessing: Feature Selection 

I now began the next stage of the preprocessing pipeline, feature selection. Feature selection is a common stage in ML/big data analysis, so I was slightly familiar with the concept going into this phase. The authors' goal was to *calculate a subset of features that exhibit high cell-to-cell variation in the dataset (i.e, they are highly expressed in some cells, and lowly expressed in others)*. This is the condition from which they selected relevant features. 

```{r, message=FALSE} 
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
```

# Preprocessing: Scaling the data
The final stage of preprocessing was scaling, which is necessary prior to applying pricipal component analysis (PCA). The authors provided the *ScaleData()* function that: 

  1. Shifts the expression of each gene, so that the mean expression across cells is 0
  2. Scales the expression of each gene, so that the variance across cells is 1
  3. This step gives equal weight in downstream analyses, so that highly-expressed genes do not dominate
  4. The results of this are stored in pbmc[["RNA"]]@scale.data
  
```{r, message=FALSE} 
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
```

# Data Analysis: PCA 

Now that preprocessing was complete, I applied PCA to reduce the data dimensionality. The below plot shows the first 15 PCs. 

```{r, message=FALSE} 
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
# Examine and visualize PCA results a few different ways
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
DimPlot(pbmc, reduction = "pca")
# display the first 15 PCs
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)
```

# Data Analysis: Dimensionality Determination 

**Seurat** clusters cells by PCA scores. Since PCs are linear combinations of variables, the authors compare PCs to *metafeatures*. The common question in PCA is how many PCs to incorporate. The authors provide 2 possible methods to answer to this question. 

```{r, message=FALSE} 
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
# Visualization 1: 
JackStrawPlot(pbmc, dims = 1:15)

# Visualization 2: 
ElbowPlot(pbmc)
```

# Data Analysis: Clustering 

```{r, message=FALSE}
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
# Look at cluster IDs of the first 5 cells
head(Idents(pbmc), 5)
```

# Data Analysis: Non-linear Dimesionality Reduction

```{r, message=FALSE}
# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
pbmc <- RunUMAP(pbmc, dims = 1:10)
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(pbmc, reduction = "umap")
```

# Data Analysis: Finding cluster biomarkers
```{r, message=FALSE}
# find all markers of cluster 2
cluster2.markers <- FindMarkers(pbmc, ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, n = 5)

# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)

# find markers for every cluster compared to all remaining cells, report only the positive ones
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

cluster0.markers <- FindMarkers(pbmc, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
VlnPlot(pbmc, features = c("MS4A1", "CD79A"))
# plot raw counts
VlnPlot(pbmc, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)
FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", 
    "CD8A"))
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(pbmc, features = top10$gene) + NoLegend()
```

# Data Analysis: Assign cell IDs to clusters
```{r, message=FALSE}
new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono", 
    "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
```
