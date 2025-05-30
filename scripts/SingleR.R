library (SingleR)

#ref = celldex::HumanPrimaryCellAtlasData()

#default for SingleR is to perform annotation of each individual cell

#Identification of monocyte subsets

#ref2 = celldex::DatabaseImmuneCellExpressionData()
View(as.data.frame(colData(ref)))
View(A_nt@meta.data)

#default for SingleR is to perform annotation of each individual cell

pbmc_counts = GetAssayData(A_nt, slot = "counts")
pbmc_counts

pred = SingleR(test = pbmc_counts, ref = ref, labels = ref$label.fine)
pred

pred2 =as.data.frame(pred)

pred2

#table = write.csv(pred, file = "labels_A.csv")
#table_2 =read.csv("labels_aqp4.csv")

pred$scores

table(pred$labels)

plotScoreHeatmap(pred)

plotDeltaDistribution(pred)

tab = table(Assigned=pred$labels, Clusters=A_nt$seurat_clusters)
tab

#tab1 = write.xlsx(tab, file = "clusters_a_subset_human_atlas.xlsx")
#tab1

pheatmap(log10(tab+10), color = colorRampPalette(c("white", "blue"))(10))

View(AQP4_nt@meta.data)

A_nt$singleR.labels = pred$labels[match(rownames(A_nt@meta.data), rownames(pred))]
A_nt$singleR.labels

DimPlot(A_nt, reduction = "umap", group.by = "singleR.labels")
View(A_nt@meta.data)

p4 = DimPlot(A_nt, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend() + ggtitle("seurat_clusters") + theme(plot.title = element_text(hjust = 0.5))
