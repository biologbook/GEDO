GEDO - Topology-based inference of gene module activity in Sjögren’s disease.

# 📘 Abstract
Inferring gene module activity is key to understanding transcriptomic dysregu-
lation in case-control studies. Most existing Gene Set Analysis methods focus on
the relative positioning of data points without considering the complex geome-
try of the data.
We introduce GEDO, a graph- and topology-based method that infers gene
module activity through a transition score quantifying the shift from healthy
controls to diseased individuals. When applied to bulk RNA-seq data from
Sjögren’s disease patients and healthy controls (PRECISESADS cohort), GEDO
was benchmarked against PCA1 and the Mean of z-scores in classification and
unsupervised clustering tasks.
GEDO outperformed both approaches in predicting disease status and yielded
more informative module matrices, enhancing subgroup separability. The bio-
logical signal captured was aligned with the knowledge and clinical features of
Sj¨ogren’s disease.
GEDO offers a robust, interpretable framework for quantifying gene modules
activity, with applications in single- and multi-omics integration. Its topology-
based design enables finer resolution of disease-related transcriptomic alterations.

# 📂 Structure of the repository
docs/TUTORIAL.html : Quarto edited file which run analysis, show code execution, and show figures.
R/functions_article.R : functions for computing module matrices with all methods exepts GEDO, and ton run comparison analyses.
R/GEDO.R : not provided in this repository. Contain functions necessary to run GEDO algorithm. This script can be shared upon request to the authors, and under certain conditions. 
scripts/GEDO_toy_example.R : script which run GEDO on random values generated directly in the script, and plot a heatmap. 
scripts/compute_GEDO.qmd : Mardown Quarto to run GEDO and get figures.

# 📄 Tutorial 
👉 [See quarto tutorial](docs/TUTORIAL.md)


