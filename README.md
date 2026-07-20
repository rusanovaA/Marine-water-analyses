# Marine Water Metagenomic Analysis

Scripts and data used for analysis and visualization of **marine water metagenomic datasets**.

---

## Repository structure

### Main scripts
- **`Main.R`**  
  Used to create a `phyloseq` object from the raw input tables and perform visualization, statistical testing, and figure generation.  

- **`Additional.R`**  
  Used to create a `phyloseq` object from the raw input tables and perform additional analyses.

---

## Input data

All source files are located in the folder **`Source_data/`**:

- `all_OTU_frequency.tsv` — OTU abundance table.  
- `all_OTUs_phylogeny.tsv` — taxonomic and phylogenetic information for OTUs.  
- `Supplementary_Table_1_St_Anna_Trough_metadata_0m_50m.tsv` — sample metadata.

These files were used to construct the **phyloseq** object and perform downstream analyses and visualizations.

---

## Network analysis (SpeSpeNet)

Files located in **`Source_data/SpeSpeNet/`** were used to generate the final **network visualizations**, based on results computed in **SpeSpeNet v1**  
([https://tbb.bio.uu.nl/SpeSpeNet](https://tbb.bio.uu.nl/SpeSpeNet)).

---

## SpeSpeNet cluster composition data

The following tables provide node and cluster composition information for all core samples, as well as separately for **0 m** and **50 m** depths:

- `Supplementary_Table_3_SpeSpeNet_Node_data_2025-11-05_0m_50m.txt`  
- `Supplementary_Table_4_SpeSpeNet_Node_data_2025-11-05_0m.txt`  
- `Supplementary_Table_5_SpeSpeNet_Node_data_2025-11-05_50m.txt`

---

## Summary

This repository contains:

- Scripts for **phyloseq object construction** and **data visualization**.  
- Input and result files for **network inference and cluster analysis** using *SpeSpeNet*.  
- Reproducible data and code used for **marine water metagenome analysis**.

---

## Reproducibility

The analysis was performed in **R (4.2.2)** using the following key packages:

- `phyloseq`  
- `ggplot2`  
- `ggrepel`  
- `dplyr`  
- `readxl`  
- `pairwiseAdonis`  
- `RColorBrewer`  
- `vegan`  *(PERMANOVA, ordination)*  
- `ggpubr`  
- `microbiome`  
- `ranacapa`  *(rare curve visualization)*  
- `MicrobeR`  *(3D PCoA plots)*  
- `umap`  
- `dendextend`  *(dendrogram visualization)*  
- `microbiomeMarker`  *(LefSe analysis)*  
- `pheatmap`  
- `gridExtra`  
- `VennDiagram`  *(Venn diagrams)*  
- `pals`  
- `tidygraph`  
- `ggraph`  
- `ggdark`  
- `igraph`

Operating system: **Windows 10**  
Last updated: **July 2026**
