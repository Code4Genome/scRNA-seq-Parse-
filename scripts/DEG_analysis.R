library(Seurat)
library(patchwork)
library(SeuratObject)
library(metap)
library(SingleR)
library(tidyr)

rm(list = ls())

SaveObject <- function(object, name){
  saveRDS(object, paste0(data_path, name, ".RDS"))
}

ReadObject <- function(name){
  readRDS(paste0(data_path, name, ".RDS"))
}

data_path = "~/Parse_analysis/AQP4_NT vs.HD/"
fig_path = "~/Parse_analysis/AQP4_NT_vs.HD/"

SaveFigure <- function(plots, name, type = "png", width, height, res){
  if(type == "png") {
    png(paste0(fig_path, name, ".", type),
        width = width, height = height, units = "in", res = 200)
  } else {
    pdf(paste0(fig_path, name, ".", type),
        width = width, height = height)
  }
  print(plots)
  dev.off()
}

nt = readRDS("a_nt.rds")
nt
hd = readRDS("h.rds")
hd


nt_vs_hd = merge(x = nt, y =hd , add.cell.ids = c("A_NT", "H"), project = "A")
nt_vs_hd@meta.data$orig.ident

nt_vs_hd@meta.data$orig.ident = factor (rep("A_vs_H", nrow(nt_vs_hd@meta.data)))
Idents(nt_vs_hd) = nt_vs_hd@meta.data$orig.ident

#cell quality control 
View(nt_vs_hd@meta.data)

nt_vs_hd[["percent.mt"]] <- PercentageFeatureSet(nt_vs_hd, pattern = "^MT-")
plot = VlnPlot(nt_vs_rtx, pt.size = 0.10,
               features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot

SaveFigure(plot, "vln_QC_subset", width = 12, height = 3)


nt_vs_hd@meta.data$sample = rownames(nt_vs_hd@meta.data)
View(nt_vs_hd@meta.data)

#split sample column

nt_vs_hd@meta.data = separate(nt_vs_hd@meta.data, col = "sample", into =c ("condition"), sep = "_")
View(nt_vs_hd@meta.data)

#calculate mitochondrial percentage 
nt_vs_hd$mitoPercent = PercentageFeatureSet(nt_vs_hd, pattern = "^MT-") 

#explore QC

plot2 = FeatureScatter(nt_vs_hd, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth(method = "lm")
plot2

#Filtering

filtered = subset(nt_vs_hd, subset = nFeature_RNA > 200 & nFeature_RNA < 6000  & nCount_RNA > 1000 & percent.mt < 10)
filtered

plot3 = VlnPlot(filtered, pt.size = 0.10,
                features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot3
SaveFigure(plot3, "vln_QC_subset_after_filtering", width = 12, height = 3)


#perform standard workflow steps to figure out if we see any batch effect

filtered = NormalizeData(object = filtered,  normalization.method = "LogNormalize", scale.factor = 1e4)
filtered = FindVariableFeatures(object = filtered, selection.method = "vst", nfeatures = 2000)
filtered = ScaleData(object = filtered, verbose = F)
filtered = RunPCA(object = filtered, npcs = 50, verbose = F)

ElbowPlot(filtered)

filtered = FindNeighbors(object = filtered, reduction = "pca" ,dims = 1:15)
filtered = FindClusters(object = filtered, resolution = 0.1)
filtered = RunUMAP(filtered, dims = 1:15)

#plot

DimPlot(filtered, reduction = "umap")
DimPlot(filtered, reduction = "umap", group.by = "condition")
DimPlot(filtered, reduction = "umap", group.by = "orig.ident")


FeaturePlot(filtered, c("HLA-DQA1", "HLA-DRB1", "HLA-DQB1", "CXCL8", "TGFBR3", "TGFBR1"), ncol = 2, pt.size = 0.1)


#perform integration to correct for batch effects

obj.list = SplitObject(filtered, split.by = "condition")
obj.list
for( i in 1: length(obj.list)){
  obj.list[[i]] = NormalizeData(object = obj.list[[i]])
  obj.list[[i]] = FindVariableFeatures(object = obj.list[[i]])
}  

#select integration features
features = SelectIntegrationFeatures(object.list = obj.list)

#find anchors to integrate the data
anchors = FindIntegrationAnchors(object.list = obj.list, anchor.features = features)

#integrate data

seurat.integrate = IntegrateData(anchorset = anchors)

#Scale Data, run PCA and UMAP and visualize integrated data
seurat.integrate = ScaleData(object = seurat.integrate)
seurat.integrate = RunPCA(object = seurat.integrate)
seurat.integrate = RunUMAP(object = seurat.integrate, dims = 1:50)

p4 = DimPlot(seurat.integrate, reduction = 'umap', group.by = "condition")
p4

p5 = DimPlot(seurat.integrate, reduction = 'umap', group.by = "singleR.labels")
p5


#Find markers between conditions 
#always check if it is RNA or not
DefaultAssay(seurat.integrate) = "RNA"
View(seurat.integrate@meta.data)

seurat.integrate@meta.data$cell.type = seurat.integrate@meta.data$singleR.labels
seurat.integrate@meta.data$cell.type = paste0(seurat.integrate$condition, "_",seurat.integrate$cell.type)
View(seurat.integrate@meta.data)

Idents(seurat.integrate) = seurat.integrate$cell.type

DimPlot(seurat.integrate, reduction = "umap", label = TRUE)

#Find markers 

pseudo_bulk = AggregateExpression(seurat.integrate, assays = "RNA",
                              return.seurat = T, group.by = "condition")

tail(Cells(pseudo_bulk))


bulk = FindMarkers(seurat.integrate, ident.1 =  "A_B_cell", ident.2 = "H_B_cell", test.use = "DESeq2")

head(bulk, n = 15)

a = write.csv(bulk, "pseudobulk.csv")

sum(bulk$p_val_adj <0.05 & bulk$avg_log2FC > 0.25, na.rm = T)
upregulated <- subset(bulk, avg_log2FC > 0.25 & p_val_adj < 0.05)
upregulated

sum(bulk$p_val_adj <0.05 & bulk$avg_log2FC < -0.25, na.rm = T)
downregulated <- subset(bulk, avg_log2FC < -0.25 & p_val_adj < 0.05)
downregulated

b = write.csv(upregulated, "pseudobulk_upregulated.csv")

c = write.csv(downregulated, "pseudobulk_downregulated.csv")
