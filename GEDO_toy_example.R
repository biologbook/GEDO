
# 1. Packages and functions-------------------------------------------------------------


folder_script_to_source = "~/these_clement/studies/R_Thèse_clement_LBAI/11_GEDO/Script/"
#folder to find scripts GEDO.R and functions_article.R.

source(file = paste0(folder_script_to_source,"GEDO.R"))
source(file = paste0(folder_script_to_source,"functions_article.R"))

detach_all_packages()
packages_list = c("data.table", "FNN", "magrittr", "tryCatchLog","rlist",
                  "igraph", "pbapply","rgl","rdist","msigdbr","msigdbdf",
                  "future.apply","dbscan", "progressr","pheatmap",
                  "gridExtra","MASS","plotly","ggplot2", "patchwork",
                  "caret","randomForest", "dplyr", "tidyr","parallel", 
                  "pROC", "ggpubr","fpc","cluster","phateR","reticulate", 
                  "devtools", "bench","RANN")
install_and_load(packages_list)



# 2. Creating theoretical data --------------------------------------------

reactome_modules <- data.table(msigdbr(species = "Homo sapiens", collection = "C7", subcollection = "IMMUNESIGDB"))
reactome_modules <- reactome_modules[, .SD, .SDcols = c("gs_name", "ensembl_gene")]
modules = unique(reactome_modules$gs_name)[1:3]
genes = unique(reactome_modules[gs_name %in% modules]$ensembl_gene)
n_genes=100
genes = genes[1:n_genes]

n=500
dt <- data.table(matrix(rnorm(n * length(genes)), nrow = n, ncol = length(genes)))
setnames(dt, genes)  
dt <- dt[, lapply(.SD, function(x) (x - min(x)) / (max(x) - min(x)))]

#Diagnosis vector
diag = c(rep("Diseased",n/2), rep("Control",n/2))

#IFN score of individuals
ifn_score = rnorm(n = n, mean = 2, sd = 1)



# 3. Computing GEDO -------------------------------------------------------

gedo_obj = gedo(data = dt, diag = diag, reference_group = "Control",
                category = "C7",subcategory = "IMMUNESIGDB",k_lof = 30,
                core_pct = 0.9, k_graph = 15, dim_reduc_method="none", num_cores=num_cores,
                distance="euclidean")


# 4. Computing heatmap of GEDO module matrix ------------------------------
print(heatmap.ifn(gedo_obj = gedo_obj, IFN_score=ifn_score))


