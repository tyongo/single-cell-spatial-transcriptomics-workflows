# single-cell-spatial-transcriptomics-workflows


***

# Single-Cell & Spatial Transcriptomics: Glioblastoma Heterogeneity & Trajectory Analysis

## Problem
This project dissects the cellular heterogeneity of Glioblastoma (GBM) tumors to resolve the spatial organization of molecular subtypes (Mesenchymal vs. Proneural) and infer the regulatory trajectories driving cell differentiation.

## Methods
This repository implements a multi-modal analysis pipeline integrating Single-Cell RNA-seq (scRNA-seq) and Spatial Transcriptomics (ST) data:

*   **scRNA-seq Processing (Seurat):** Dimensionality reduction (PCA/UMAP) and unsupervised clustering to identify distinct cell populations.
*   **Biomarker Discovery:** Differential expression analysis to characterize cluster-specific markers (e.g., *SPARCL1*, *FABP7*) and functional association mapping using the STRING database.
*   **Spatial Domain Identification (DR-SC):** Implementation of a Dimension-Reduction Spatial-Clustering (DR-SC) framework using a Hidden Markov Random Field (HMRF) model to enforce spatial smoothness and segment tissue domains (e.g., cortical layers or tumor boundaries).
*   **Trajectory Inference & Gene Regulation (SINCERA/Slingshot):** Construction of pseudotime lineages to model cell fate decisions and inference of Transcriptional Regulatory Networks (TRN) to identify "driving force" transcription factors.

## Results
*   **Cluster Annotation Table:** `results/tables/GBM_cluster_markers.csv` (Contains log2FC and p-values for top markers including *SPARCL1* and *FABP7*).
*   **Dimensionality Reduction:** `results/figures/UMAP_clusters_feature_plots.png` (Visualizing expression of *TTYH1*, *SLN*, and *GDF1*).
*   **Spatial Segmentation:** `results/figures/Spatial_heatmap_DRSC.png` (Spatially resolved tumor subtypes).
*   **Regulatory Network:** `results/figures/TRN_driving_force.pdf` (Directed graph of TF-target interactions).

## Reproducibility

```bash
# Clone the repository
git clone https://github.com/terver-yongo/single-cell-spatial-transcriptomics-workflows.git
cd single-cell-spatial-transcriptomics-workflows

# Create the environment from the provided YAML file
conda env create -f environment.yml

# Activate the environment
conda activate sc_spatial_env

# Run the Seurat analysis pipeline (R)
Rscript notebooks/01_GBM_scRNA_Seurat_Pipeline.Rmd

# Run the Spatial Clustering pipeline (DR-SC)
Rscript notebooks/02_Spatial_Domain_ID_DR-SC.R

# Run Trajectory Inference
Rscript notebooks/03_Trajectory_TRN_Inference.Rmd
```

## Skills Demonstrated

### Spatial & Single-Cell Analysis
Demonstrated proficiency in **Seurat** and **Scanpy** for high-dimensional data processing. Applied advanced algorithms like **DR-SC** (joint dimension reduction and spatial clustering) to overcome the limitations of standard PCA in spatially resolved data, successfully identifying biologically relevant tissue domains.

### Translational Omics & Biomarker Discovery
Executed rigorous differential expression workflows to identify diagnostic markers (*SPARCL1*, *FABP7*) with high log-fold changes. Integrated external databases (**STRING**) to validate functional protein-protein interactions (e.g., SPARCL1-MYOC association), translating computational clusters into biological mechanisms relevant to Glioblastoma subtypes.

### Systems Biology & Network Inference
Utilized **SINCERA** and **Slingshot** to move beyond static clustering into dynamic modeling. Constructed **Transcriptional Regulatory Networks (TRNs)** to identify key regulatory "hubs" and inferred pseudotime trajectories to model the developmental lineage of tumor cells.

## Next Improvements
*   **Multi-sample Integration:** Implement **Harmony** or **Canonical Correlation Analysis (CCA)** to integrate datasets across multiple patients and correct for batch effects.
*   **Deep Learning Integration:** Incorporate **stAI** (Spatial Transcriptomics AI) to impute missing genes and enhance cell-type annotation resolution by leveraging reference scRNA-seq embeddings.
*   **Interactive Visualization:** Develop a **Shiny** app or **Cellxgene** instance to allow lay-users to explore the spatial domains and gene expression patterns interactively.
