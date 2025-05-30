library(Seurat)
library(harmony)
library(DoubletFinder)
library(tidyverse)
library(ggplot2)

S1 = readRDS("S1_new.rds")
S2 = readRDS("S2_new.rds")
S3 = readRDS("S3_new.rds")

S4 = readRDS("S1_2new.rds")
S5 = readRDS("S2_2new.rds")
S6 = readRDS("S3_2new.rds")

#merge data

merged_seurat = merge(S1, y = c(S2, S3, S4, S5, S6), add.cell.ids = c("s1", "s2", "s3", "s4", "s5", "s6"), project = "NMO")

merged_seurat@meta.data$orig.ident

merged_seurat@meta.data$orig.ident = factor (rep("NMO", nrow(merged_seurat@meta.data)))
Idents(merged_seurat) = merged_seurat@meta.data$orig.ident

str(merged_seurat)

#QC and filtering

merged_seurat[["percent.mt"]] <- PercentageFeatureSet(merged_seurat, pattern = "^MT-")
plot = VlnPlot(merged_seurat, pt.size = 0.10,
               features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot

plot2 = FeatureScatter(merged_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth(method = "lm")
plot2

#Filtering 

filtered2 = subset(merged_seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & nCount_RNA < 25000 & percent.mt < 15)
filtered2
plot3 = VlnPlot(filtered2, pt.size = 0.10,
                features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot3

#performs standard workflow

debug(sum)

#standard workflow steps
filtered2 = NormalizeData(filtered2, normalization.method = "LogNormalize", scale.factor = 10000)
filtered2 = FindVariableFeatures(filtered2, selection.method = "vst", nfeatures = 2000)
filtered2 = ScaleData(filtered2)
filtered2 = RunPCA(filtered2)
ElbowPlot(filtered2, ndims = 50)

filtered2 = FindNeighbors(filtered2, dims =1:30)
filtered2 = FindClusters(filtered2)
filtered2 = RunUMAP(filtered2, dims = 1:30)

before = DimPlot(filtered2, reduction = "umap", group.by ="sample")

### pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.list_data <- paramSweep(filtered2, PCs = 1:30, sct = FALSE)
sweep.stats_data <- summarizeSweep(sweep.res.list_data, GT = FALSE)
bcmvn_data <- find.pK(sweep.stats_data)

ggplot(bcmvn_data, aes(pK, BCmetric, group = 1)) +
  geom_point() +
  geom_line()

#select the pK value

pK = bcmvn_data %>%
  filter(BCmetric == max(BCmetric)) %>%
  select(pK)

pK = as.numeric(as.character(pK[[1]]))

annotations = filtered2@meta.data$seurat_clusters

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.03*nrow(filtered2@meta.data))  ## Assuming 3% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

filtered2 <- doubletFinder(filtered2, 
                           PCs = 1:30, 
                           pN = 0.25, 
                           pK = pK, 
                           nExp = nExp_poi.adj, 
                           reuse.pANN = "FALSE", 
                           sct = FALSE)

View(filtered2@meta.data$DF.classifications_0.25_0.01_203)
#visualize
DimPlot(filtered2, reduction = "umap", group.by ="DF.classifications_0.25_0.01_203")

#number of singlets and doublets
table(filtered2@meta.data$DF.classifications_0.25_0.01_203)

filtered2 = subset(filtered2, subset = DF.classifications_0.25_0.01_203 == "Singlet")

table(filtered2@meta.data$DF.classifications_0.25_0.01_203)

#plot

#Run Harmony

harmony <- filtered2 %>%
  RunHarmony(group.by.vars = "sample", plot_convergence = TRUE)

harmony@reductions

harmony.embeddings = Embeddings(harmony, "harmony")
harmony.embeddings[1:10, 1:10]

p1 <- DimPlot(object = harmony, reduction = "harmony", pt.size = .1, group.by = "sample") + NoLegend()

harmony = harmony %>%
  RunUMAP(reduction = "harmony", dims = 1:30) %>%
  FindNeighbors(reduction = "harmony", dims = 1:30) %>%
  FindClusters(resolution = 0.5)

after = DimPlot(harmony, reduction = "umap", group.by = "sample")

before|after

saveRDS(harmony, file = "NMO_6samples_harmony.rds")

#-------------------------------------------------------------

harmony= BuildClusterTree(harmony, reorder = TRUE, reorder.numeric = TRUE)

markers = FindAllMarkers(harmony, min.pct = 0.25, logfc.threshold = 0.25)
markers
markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)


top =  markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
top
to_plot = unique(top$gene)
to_plot

plot4 = DotPlot(harmony, features = to_plot, group.by = "tree.ident") + coord_flip()
plot4

DoHeatmap(harmony, features = top$gene, group.bar.height = 0.01, size = 3) + NoLegend()
View(harmony@meta.data)

ref = celldex::MonacoImmuneData()
View(as.data.frame(colData(ref)))
View(combined@meta.data)

harmony@assays$RNA@data
