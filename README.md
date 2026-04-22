# 📘 GEDO
GEDO: Topology-based inference of gene module activity in Sjögren’s disease.
doi : 

## Abstract
Inferring gene module activity is key to understanding transcriptomic dysregulation in case-control studies. Most existing Gene Set Analysis methods focus on differences between groups without considering the complex geometry of the data.

We introduce GEDO, a graph and topology-based method that infers gene module activity through a transition score, quantifying the shift from healthy controls to diseased individuals. When applied to bulk RNA-seq data from Sjögren’s disease patients and healthy controls (PRECISESADS cohort), and a breast cancer dataset from The Cancer Genome Atlas, GEDO was benchmarked against PCA, the mean of z-scores, ssGSEA, and GSVA in classification, unsupervised clustering tasks and robustness against noise and bias.

On the PRECISESADS cohort, GEDO outperformed the other approaches in predicting disease status, interferon signature, enhancing subgroup separability, and robustness against noise and bias. The biological signal captured was aligned with the knowledge and clinical features of Sjögren’s disease. In the breast cancer dataset, GEDO's embeddings better represent PAM50 molecular subtypes. 

Its supervised and topology-based design enables finer resolution of disease-related transcriptomic alterations. GEDO offers a robust, interpretable framework for quantifying gene modules' activity, with applications in single- and potentially in multi-omics integration.




# 📘 mGEDO

Alpha-Interferons and CD40: complementary pathways revealed by mGEDO in distinct transcripto-methylomic profiles of patients with Sjögren’s disease.
doi :

## Abstract : 

Objectives: Sjögren’s disease (SjD) is clinically and molecularly heterogeneous, which limits robust patient stratification and may contribute to variable responses to targeted therapies. The objective of this study is to quantify alpha-Interferon, CD40 and Interleukin 2 pathways in different SjD subgroups.

Methods: We developed mGEDO, a multi-omic extension of GEDO, to integrate whole-blood RNA-seq and DNA methylation data at the gene-set level using BloodGen3 modules. We applied mGEDO to 702 PRECISESADS participants to get two-dimensional embeddings. We then designed weighted alpha-Interferon-, CD40-, and IL2-related signatures using an automated gene selection method, and computed patient-level scores using z-score summaries and customized regularized canonical correlation analysis.

Results: mGEDO captured a continuous SjD–control molecular landscape and revealed substantial within-group heterogeneity beyond discrete clustering. Alpha-Interferon, CD40, and IL-2 signatures showed distinct but partially overlapping patterns. Alpha-Interferon remained the dominant axis associated with severe biological activity, while CD40 and IL-2 provided complementary information, particularly in IFN-negative subsets. Stratification by combined IFN/CD40 status identified biologically meaningful patient classes, with double-positive patients displaying the strongest immune dysregulation.

Conclusions: mGEDO provides an interpretable framework for integrating transcriptomic and methylomic data and for mapping SjD as a molecular continuum. The complementarity between IFN and CD40/IL2 signatures supports partially uncoupled innate and adaptive immune programs and may improve patient stratification for precision therapeutic targeting.




# 📂 Repository structure

```
GEDO/
├── README.md
├── LICENSE
├── GEDO.Rproj
│
├── R/
│   ├── GEDO/
│   │   ├── functions_article.R
│   │   │   Functions for computing module matrices with all methods except GEDO,
│   │   │   and to run comparison analyses.
│   │   └── GEDO.R
│   │       Not provided in this repository. Contains functions necessary to run the GEDO algorithm. This script can be shared upon request to the authors, under certain conditions. It must be placed in the R/GEDO/ folder.
│   │
│   └── mGEDO/
│       ├── mGEDO_functions.R
│           Functions implementing the mGEDO method (extension of GEDO),
│           including modified algorithms and additional analysis utilities.
│       └── mGEDO.R
│           Not provided in this repository. Contains functions necessary to run the mGEDO algorithm. This script can be shared upon request to the authors, under certain conditions. It must be placed in the R/mGEDO/ folder.
│
├── article_GEDO/
│   ├── README.md
│   │
│   ├── docs/
│   │   └── TUTORIAL.md
│   │       Quarto edited file which runs analysis, shows code execution, and displays figures.
│   │
│   └── scripts/
│       ├── compute_GEDO.qmd / compute_GEDO.rmarkdown
│       │   Quarto Markdown file to run GEDO and generate figures.
│       ├── GEDO_toy_example.R
│       │   Script which runs GEDO on randomly generated values and plots a heatmap.
│       ├── 1_grid_search.R
│       │   Script to run GEDO on each configuration tested for sensitivity analysis.
│       ├── 2_evaluation_grid.R
│       │   Script to compute classification metrics for each configuration.
│       ├── 3_cancer_subtypes.R
│       │   Script to run all analyses on the TCGA-BRCA dataset.
│       ├── 4_functions_reviewing_TCGA.R
│       │   Functions required for TCGA-BRCA dataset analysis.
│       ├── 5_data_points_needed.R
│       │   Script to analyze GEDO performance depending on the number of observations.
│       ├── 6_clustering.R
│       │   Script for more detailed clustering analysis than in compute_GEDO.qmd.
│       └── 7_classification.R
│           Script for more detailed classification analysis than in compute_GEDO.qmd.
│
├── article_mGEDO/
│   ├── README.md
│   │
│   └── scripts/
│       └── mGEDO_compute.qmd
│           Quarto workflow to run the mGEDO method and generate associated results.
```



# 📄 Tutorial 
## GEDO
👉 [See quarto tutorial](article_GEDO/docs/TUTORIAL.md)
[▶️ See explicative video](https://github.com/biologbook/GEDO/releases/tag/uploading_GEDO_manim_animation/GEDO_manim_animation.mp4)

