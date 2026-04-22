GEDO - Topology-based inference of gene module activity in Sjögren’s disease.

doi : 

# 📘 Abstract
Inferring gene module activity is key to understanding transcriptomic dysregulation in case-control studies. Most existing Gene Set Analysis methods focus on differences between groups without considering the complex geometry of the data.

We introduce GEDO, a graph and topology-based method that infers gene module activity through a transition score, quantifying the shift from healthy controls to diseased individuals. When applied to bulk RNA-seq data from Sjögren’s disease patients and healthy controls (PRECISESADS cohort), and a breast cancer dataset from The Cancer Genome Atlas, GEDO was benchmarked against PCA, the mean of z-scores, ssGSEA, and GSVA in classification, unsupervised clustering tasks and robustness against noise and bias.

On the PRECISESADS cohort, GEDO outperformed the other approaches in predicting disease status, interferon signature, enhancing subgroup separability, and robustness against noise and bias. The biological signal captured was aligned with the knowledge and clinical features of Sjögren’s disease. In the breast cancer dataset, GEDO's embeddings better represent PAM50 molecular subtypes. 

Its supervised and topology-based design enables finer resolution of disease-related transcriptomic alterations. GEDO offers a robust, interpretable framework for quantifying gene modules' activity, with applications in single- and potentially in multi-omics integration.


# 📂 Repository structure
├── docs/
│   └── TUTORIAL.html
│       └── Quarto edited file which runs analysis, shows code execution, and displays figures.
│
├── R/
│   ├── functions_article.R
│   │   └── Functions for computing module matrices with all methods except GEDO, and to run comparison analyses.
│   └── GEDO.R
│       └── Not provided in this repository. Contains functions necessary to run the GEDO algorithm.
│           This script can be shared upon request to the authors, under certain conditions.
│           It must be placed in the R/ folder.
│
├── scripts/
│   ├── GEDO_toy_example.R
│   │   └── Script which runs GEDO on randomly generated values and plots a heatmap.
│   ├── compute_GEDO.qmd
│   │   └── Quarto Markdown file to run GEDO and generate figures.
│   ├── 1_grid_search.R
│   │   └── Script to run GEDO on each configuration tested for sensitivity analysis.
│   ├── 2_evaluation_grid.R
│   │   └── Script to compute classification metrics for each configuration.
│   ├── 3_cancer_subtypes.R
│   │   └── Script to run all analyses on the TCGA-BRCA dataset.
│   ├── 4_functions_reviewing_TCGA.R
│   │   └── Functions required for TCGA-BRCA dataset analysis.
│   ├── 5_data_points_needed.R
│   │   └── Script to analyze GEDO performance depending on the number of observations.
│   ├── 6_clustering.R
│   │   └── Script for more detailed clustering analysis than in compute_GEDO.qmd.
│   └── 7_classification.R
│       └── Script for more detailed classification analysis than in compute_GEDO.qmd.




# 📄 Tutorial 
👉 [See quarto tutorial](docs/TUTORIAL.md)
[▶️ See explicative video](https://github.com/user/repo/releases/download/v1.0/GEDO_manim_animation.mp4)

