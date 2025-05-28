library(Seurat)
library(sctransform)
library(tidyverse)

S1 = readRDS("S1_new.rds")
S2 = readRDS("S2_new.rds")
S3 = readRDS("S3_new.rds")

S4 = readRDS("S1_2new.rds")
S5 = readRDS("S2_2new.rds")
S6 = readRDS("S3_2new.rds")


S1[["percent.mt"]]  <- PercentageFeatureSet(S1, pattern = "^MT-")
S2[["percent.mt"]]  <- PercentageFeatureSet(S2, pattern = "^MT-")
S3[["percent.mt"]]  <- PercentageFeatureSet(S3, pattern = "^MT-")

S4[["percent.mt"]]  <- PercentageFeatureSet(S4, pattern = "^MT-")
S5[["percent.mt"]]  <- PercentageFeatureSet(S5, pattern = "^MT-")
S6[["percent.mt"]]  <- PercentageFeatureSet(S6, pattern = "^MT-")

VlnPlot(S1, features = c("nFeature_RNA","nCount_RNA","percent.mt"), ncol = 3)

VlnPlot(S2, features = c("nFeature_RNA","nCount_RNA","percent.mt"), ncol = 3)

VlnPlot(S3, features = c("nFeature_RNA","nCount_RNA","percent.mt"), ncol = 3)

VlnPlot(S4, features = c("nFeature_RNA","nCount_RNA","percent.mt"), ncol = 3)

VlnPlot(S5, features = c("nFeature_RNA","nCount_RNA","percent.mt"), ncol = 3)

VlnPlot(S6, features = c("nFeature_RNA","nCount_RNA","percent.mt"), ncol = 3)


S1 <- subset(S1, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & nCount_RNA < 20000 & percent.mt < 15)
S2 <- subset(S2, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & nCount_RNA < 20000 & percent.mt < 15)
S3 <- subset(S3, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & nCount_RNA < 20000 & percent.mt < 15)

S4 <- subset(S4, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & nCount_RNA < 20000 & percent.mt < 15)
S5 <- subset(S5, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & nCount_RNA < 20000 & percent.mt < 15)
S6 <- subset(S6, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & nCount_RNA < 20000 & percent.mt < 15)

#-----------------------------------------------------------------------------------------------------------

obj1 <- SCTransform(S1, vst.flavor = "v2", vars.to.regress = "percent.mt")
obj2 <- SCTransform(S2, vst.flavor = "v2", vars.to.regress = "percent.mt")
obj3 <- SCTransform(S3, vst.flavor = "v2", vars.to.regress = "percent.mt")

obj4 <- SCTransform(S4, vst.flavor = "v2", vars.to.regress = "percent.mt")
obj5 <- SCTransform(S5, vst.flavor = "v2", vars.to.regress = "percent.mt")
obj6 <- SCTransform(S6, vst.flavor = "v2", vars.to.regress = "percent.mt")

#-----------------------------------------------------------------------------------

sample.list <- list(obj1, obj2, obj3, obj4, obj5, obj6)

features <- SelectIntegrationFeatures(object.list = sample.list, nfeatures = 3000)

future::plan("sequential")
sample.list <- PrepSCTIntegration(object.list = sample.list, anchor.features = features, verbose = TRUE)

future::plan("multisession")
anchors <- FindIntegrationAnchors(object.list = sample.list, normalization.method = "SCT", anchor.features = features)
combined<- IntegrateData(anchorset = anchors, normalization.method = "SCT")

combined@assays$SCT@data

combined<- RunPCA(combined)
ElbowPlot(combined, ndims = 50)

combined<- FindNeighbors(combined, dims = 1:20)
combined<- FindClusters(combined, resolution = 0.5)

combined<- RunUMAP(combined, dims = 1:20)

after = DimPlot(combined, reduction = "umap", group.by = "sample")

before|after

saveRDS(combined, "combined_6samples.rds")

combined = readRDS("combined_6samples.rds")

