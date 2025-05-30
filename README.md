# 🧬 Single-Cell RNA-seq Analysis (Parse Biosciences) 

![grafik](https://github.com/user-attachments/assets/febcee89-562f-4ea5-a035-9396326fa867)

Parse Biosciences technology offers several key advantages for single-cell sequencing, including increased scalability, enhanced flexibility with sample fixation, and simplified workflows. 
It allows researchers to fix and store samples, then run them together in a single workflow, enabling large-scale experiments without the limitations of traditional microfluidics-based approaches. 
Also, this technology also provides combinatorial barcoding to process up to 96 samples or one million cells in a single experiment and the compatibility with both short- and long-read applications.
It offers high-quality data, long-term sample storage, and reduced time-course experiments

![grafik](https://github.com/user-attachments/assets/3e008beb-b84f-43e7-a6a1-71a8b4d1a861)


More info: www.parsebiosciences.com

## 🧬 Workflow Overview

1. **Preprocessing, Demultiplexing, UMI counting and Expression matrix Generation**
   - Input: sublibrary 
   - Tool: `parse conda environment`
   - Output: count_matrix.mtx, cell_metadata.csv, all_genes.csv

2. **Quality Control**
   - Tools: Seurat / Scanpy
   - Metrics: nFeature_RNA, nCount_RNA, mitochondrial content

3. **Normalization & Integration**
   - SCTransform, log-normalization
   - Batch effect correction (Harmony / Seurat integration)

4. **Dimensionality Reduction & Clustering**
   - PCA, UMAP

5. **Differential Expression & Marker Identification**
   - FindMarkers (Seurat), rank_genes_groups (Scanpy)

6. **Cell Type Annotation**
   - Manual curation or reference-based annotation (SingleR, CellTypist)


## 📁 Repository Structure

├── data/ # Raw and processed data (sublibrary analysis using conda environment)       
├── scripts/ # Data analysis scripts in R (Seurat, Harmony, scTransform, DEG)        
├── README.md # Project overview and instructions         
└── LICENSE

## 🧪 Project Overview

**Assay:** Whole transcriptome single-cell RNA-seq    
**Sample Type:** Human PBMCs
