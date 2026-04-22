
# 1.Packages -------------------------------------------------------------


folder_script_to_source = "/home/clem/GEDO/R/"
source(file = paste0(folder_script_to_source,"GEDO.R"))
source(file = paste0(folder_script_to_source,"functions_article.R"))
source("/home/clem/GEDO/scripts/4_functions_reviewing_TCGA.R")


detach_all_packages()
packages_list = c("data.table", "FNN", "magrittr", "tryCatchLog","rlist",
                  "igraph",
                  "pbapply","rgl","rdist","msigdbr","msigdbdf",
                  "future.apply","dbscan", "progressr","pheatmap",
                  "gridExtra","MASS","plotly","ggplot2", "patchwork",
                  "caret","randomForest", "dplyr", "tidyr","parallel",
                  "pROC", "ggpubr","fpc","cluster","phateR","reticulate",
                  "devtools", "bench","RANN",
                  # "GSVA",
                  "BiocParallel", "visdat",
                  "magick", "ggbreak","cowplot","grid", "SummarizedExperiment",
                  "TCGAbiolinks","dplyr","genefu","DESeq2")

# install_and_load(packages_list)
sapply(packages_list, require, character.only = TRUE)

# library = "immunsigdb"
# library = "hallmark"
# library="reactome"
# library="kegg"
library="oncogenic"

cat(library,"\n")

if(library=="immunsigdb"){
category = "C7"
subcategory = "IMMUNESIGDB"
}

if(library=="hallmark"){
  category = "H"
  subcategory = NULL
}

if(library=="reactome"){
  category = "C2"
  subcategory = "CP:REACTOME"
}
if(library=="kegg"){
  category = "C2"
  subcategory = "CP:KEGG_MEDICUS"
}

if(library=="oncogenic"){
  category = "C6"
  subcategory = NULL
}


folder_for_res = "/home/clem/GEDO/results/oncogenic/"

num_cores=5



# 2. RNAseq data download -----------------------------------------------------------------


get_tcga_rna_seq_data = function(){
query_rna <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

GDCdownload(query_rna)

brca_rna <- GDCprepare(query_rna)


#--------Keeping only Primary Tumor per patient, aliquot A 
coldata <- colData(brca_rna) %>% as.data.frame()

coldata <- data.frame(
  barcode = colnames(brca_rna),
  shortLetterCode = brca_rna$shortLetterCode,
  sample_type = brca_rna$sample_type,
  stringsAsFactors = FALSE
)

coldata$patient_id <- substr(coldata$barcode, 1, 12)

coldata_one <- coldata %>%
  dplyr::filter(shortLetterCode == "TP") %>%
  dplyr::arrange(barcode) %>%
  dplyr::group_by(patient_id) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup()



brca_one_sample <- brca_rna[, coldata_one$barcode]

#gene counts : 
counts <- assay(brca_one_sample, "unstranded")
dim(counts)
anyNA(counts)   # doit être FALSE




# 3. VST Transformation ---------------------------------------------------

library(DESeq2)

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = colData(brca_one_sample),
  design = ~ 1
)

dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)


vst_mat <- assay(vst(dds, blind = TRUE))
genes = rownames(vst_mat)
ids = colnames(vst_mat)

vst_mat_t = t(vst_mat)
vst_dt = data.table(vst_mat_t)

ids = substr(ids, 1, 12)
vst_dt$ID = ids
vst_dt <- vst_dt[, c("ID", setdiff(names(vst_dt), "ID")), with = FALSE]

#Standardisation 
num_cols <- setdiff(names(vst_dt), "ID")

# vst_dt[, (num_cols) := lapply(.SD, function(x) {
#   (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
# }), .SDcols = num_cols]
vst_dt[, (num_cols) := lapply(.SD, function(x) {
  rng <- range(x, na.rm = TRUE)
  if (diff(rng) == 0) {
    return(rep(0, length(x)))
  }
  (x - rng[1]) / diff(rng)
}), .SDcols = num_cols]


return(vst_dt)
}


if(file.exists("/home/clem/GEDO/data/TCGA_BRCA_RNAseq_STARvst.rds")){
  rna_seq_data = readRDS("/home/clem/GEDO/data/TCGA_BRCA_RNAseq_STARvst.rds")
}else{
  rna_seq_data = get_tcga_rna_seq_data()
  saveRDS(rna_seq_data,"/home/clem/GEDO/data/TCGA_BRCA_RNAseq_STARvst.rds")
}

#Remove version of ensembl codes in TCGA data
colnames(rna_seq_data)=sub("\\..*$", "", colnames(rna_seq_data))
  




# 4. Loading metadata -----------------------------------------------------

library(TCGAbiolinks)
pam50_df <- TCGAquery_subtype(tumor = "BRCA")[ , c("patient", "BRCA_Subtype_PAM50")]
pam50 = data.table(pam50_df)
colnames(pam50)[1] = "ID"




# 5. Visualisation of global RNAseq data ----------------------------------


ids = rna_seq_data$ID
rna_seq_data = rna_seq_data[, -1]

coldata <- data.table(readRDS("/home/clem/GEDO/data/coldata.rds"))
setnames(x = coldata,old="patient",new="ID")




# 6. computing module matrices --------------------------------------------

#Diag vector : 

coldata <- data.table(readRDS("/home/clem/GEDO/data/coldata.rds"))
#Comme groupe ctrl, on va prendre les definition=="Solid Tissue Normal" (sein sain). n=113.
#Et on va les comparer aux cancers solides (toujours du sein). et on va voir si on retrouve les sous-groupes
setnames(x = coldata, old="patient",new="ID")
tab = data.table(ID = ids)
tab[coldata, diag:=i.shortLetterCode, on="ID"]
tab[diag=="TM",diag:="TP"]
diag=tab$diag


# #to test 
# data=copy(rna_seq_data)
# reference_group = "NT"
# category = "C7"
# subcategory = "IMMUNESIGDB"
# k_lof = 30
# core_pct = 0.9
# k_graph = 15
# dim_reduc_method="none"
# distance="euclidean"
# ncomp=NULL
# dim_reduc_dist_method=NULL


#GEDO : 
cat("Run GEDO \n")
if(file.exists(paste0(folder_for_res, "gedo_obj.rds"))){
  gedo_obj = readRDS(paste0(folder_for_res, "gedo_obj.rds"))
}else{
  gedo_obj = gedo(data = rna_seq_data, diag = diag, reference_group = "NT",
                  category = category,subcategory = subcategory,k_lof = 30,
                  core_pct = 0.9, k_graph = 15, dim_reduc_method="none", num_cores=num_cores,
                  distance="euclidean")
  
  saveRDS(gedo_obj, file=paste0(folder_for_res, "gedo_obj.rds"))
}



#GEDOcorr
cat("Run GEDOcorr \n")

if(file.exists(paste0(folder_for_res, "gedo_corr_obj.rds"))){
  gedo_corr_obj = readRDS(paste0(folder_for_res, "gedo_corr_obj.rds"))
}else{
  gedo_corr_obj = gedo(data = rna_seq_data, diag = diag, reference_group = "NT",
                       category = category,subcategory = subcategory,k_lof = 30,
                       core_pct = 0.9, k_graph = 15, dim_reduc_method="none", num_cores=num_cores,
                       distance="correlation")
  
  saveRDS(gedo_corr_obj, file=paste0(folder_for_res, "gedo_corr_obj.rds"))
}



#UMAP_GEDO
cat("UMAP GEDO \n")

if(file.exists(paste0(folder_for_res, "umap_gedo_obj.rds"))){
  umap_gedo_obj = readRDS(paste0(folder_for_res, "umap_gedo_obj.rds"))
}else{
  umap_gedo_obj = gedo(data = rna_seq_data, diag = diag, reference_group = "NT",
                       category = category,subcategory = subcategory,k_lof = 30,
                       core_pct = 0.9, k_graph = 15, dim_reduc_method="umap", ncomp = 10,
                       dim_reduc_dist_method = "correlation", num_cores=num_cores,
                       distance="euclidean")
  
  saveRDS(umap_gedo_obj, file=paste0(folder_for_res, "umap_gedo_obj.rds"))
}



#PCA1
cat("Run PCA1 \n")

if(file.exists(paste0(folder_for_res, "pca1_obj.rds"))){
  pca1_module_matrix = readRDS(paste0(folder_for_res, "pca1_obj.rds"))
}else{
  pca1_module_matrix = compute_module_matrix(method="pca1", data=rna_seq_data, 
                                             diag=diag, 
                                             reference_group="NT",
                                             category=category,
                                             subcategory = subcategory,
                                             num_cores=num_cores)
  
  saveRDS(pca1_module_matrix, file=paste0(folder_for_res, "pca1_obj.rds"))
}



#Mean Z-scores
cat("Run MZS \n")

if(file.exists(paste0(folder_for_res, "mzscore_obj.rds"))){
  mean_z_score_module_matrix=readRDS(paste0(folder_for_res, "mzscore_obj.rds"))
}else{
  mean_z_score_module_matrix = compute_module_matrix(method="mean_z_score", data=rna_seq_data, 
                                                     diag=diag, 
                                                     reference_group="NT",
                                                     category=category,
                                                     subcategory = subcategory,
                                                     num_cores=num_cores)
  
  saveRDS(mean_z_score_module_matrix, file=paste0(folder_for_res,"mzscore_obj.rds"))
}




#ssGSEA : 
cat("Run ssGSEA \n")

expr_mat <- t(as.matrix(rna_seq_data))

gene_sets_gsea=data.table(msigdbr(species = "Homo sapiens", collection = category, subcollection = subcategory))
geneSets <- split(gene_sets_gsea$ensembl_gene, f = gene_sets_gsea$gs_name)

#Filtrage des genesets avec les gènes qu'on a dans les données
geneSets_filtered <- lapply(geneSets, intersect, rownames(expr_mat))



if(file.exists(paste0(folder_for_res, "ssgsea_obj.rds"))){
  ssgsea_obj=readRDS(paste0(folder_for_res, "ssgsea_obj.rds"))
}else{
  #Paramétrage ssGSEA
  param_ssgsea <- ssgseaParam(
    exprData    = expr_mat,
    geneSets    = geneSets_filtered,
    minSize     = 3,
    normalize=T
  )
  
  scores_ssgsea  <- gsva(param_ssgsea, BPPARAM=MulticoreParam(workers = num_cores,progressbar = T))
  scores_ssgsea=data.table(t(scores_ssgsea))
  scores_ssgsea[, (names(scores_ssgsea)) := lapply(.SD, function(x) (x - min(x, na.rm=TRUE)) / (max(x, na.rm=TRUE) - min(x, na.rm=TRUE)))]
  
  
  ssgsea_obj=list(
    module_matrix = scores_ssgsea,
    config=list(param_ssgsea=param_ssgsea, reference_group="NT"),
    diag=as.factor(diag)
  )
  
  saveRDS(ssgsea_obj, file = paste0(folder_for_res, "ssgsea_obj.rds"))
}



#GSVA
cat("Run GSVA \n")

if(file.exists(paste0(folder_for_res, "gsva_obj.rds"))){
  gsva_obj=readRDS(paste0(folder_for_res, "gsva_obj.rds"))
}else{
  
  param_gsva <- gsvaParam(
    exprData               = expr_mat,            
    geneSets               = geneSets_filtered,    
    minSize                = 3,                   
    kcdf                   = "Gaussian",           
    tau                     = 1,                  
    maxDiff                = TRUE,                
    absRanking             = FALSE
  )
  
  scores_gsva  <- gsva(param_gsva, BPPARAM=MulticoreParam(workers = num_cores,progressbar = T))
  
  scores_gsva=data.table(t(scores_gsva))
  scores_gsva[, (names(scores_gsva)) := lapply(.SD, function(x) (x - min(x, na.rm=TRUE)) / (max(x, na.rm=TRUE) - min(x, na.rm=TRUE)))]
  
  gsva_obj=list(
    module_matrix = scores_gsva,
    config=list(param_gsva=param_gsva, reference_group="NT"),
    diag=as.factor(diag)
  )
  
  saveRDS(gsva_obj, file = paste0(folder_for_res, "gsva_obj.rds"))
}


cat("Computing GEDO obj : Finished ! \n")



# 7. Plot PHATE ~ subgroups -----------------------------------------------------------


module_matrix_list = list(GEDO=gedo_obj$module_matrix,
                   GEDOcorr = gedo_corr_obj$module_matrix,
                   UMAP_GEDO=umap_gedo_obj$module_matrix,
                   PCA1=pca1_module_matrix$module_matrix,
                   MZS=mean_z_score_module_matrix$module_matrix,
                   ssGSEA = ssgsea_obj$module_matrix,
                   GSVA = gsva_obj$module_matrix)


library(data.table)
library(ggplot2)
library(patchwork)
library(uwot)
library(phateR)






cat("UMAP & PHATE comparisons... \n")


if(file.exists(paste0(folder_for_res, "umap_phate_methods_comparison.pdf"))==F){

plot_list <- Map(
  function(mat, nm) {
    message("Processing method: ", nm)
    make_umap_phate_plots_tcga(
      rna_seq_data = mat,
      ids = ids,
      pam50 = pam50,
      coldata = coldata,
      method_name = nm,
      metric = "correlation"
    )
  },
  module_matrix_list,
  names(module_matrix_list)
)


plot_list_flat <- unlist(plot_list, recursive = FALSE)
final_plot <- wrap_plots(plot_list_flat, ncol = 4)

final_plot +
  plot_annotation(
    tag_levels = "A"
  )


# final_plot <- wrap_plots(plot_list, ncol = 1)
# 
# final_plot <- Reduce(`/`, plot_list)

final_plot

ggsave(plot=final_plot, filename=paste0(folder_for_res, "umap_phate_methods_comparison.pdf"), height = 30, width=30)
}




# 8. Classification  ------------------------------------------------------



library(data.table)
library(randomForest)
library(pROC)
library(ggplot2)




#Ajouter la colonne ID 
module_matrix_list <- lapply(
  module_matrix_list,
  function(dt) {
    dt[, ID := ids]        # ajoute la colonne
    setcolorder(dt, "ID")  # la met en première position
    dt
  }
)
saveRDS(module_matrix_list, file = paste0(folder_for_res,"module_matrix_list.rds"))
# cat("Compute umap of each module matrix \n")
# module_matrix_list_umap <- lapply(
#   module_matrix_list,
#   function(dt) {
#     set.seed(123)
#     umap = data.table(uwot::umap(X = dt,n_neighbors = 15,n_components = 2,metric = "correlation"))
#     colnames(umap) = c("UMAP1","UMAP2")
#     umap[, ID := ids]        # ajoute la colonne
#     setcolorder(umap, "ID")  # la met en première position
#     return(umap)
#   }
# )



cat("Classification... \n")


if(file.exists(paste0(folder_for_res,"KNN_results.rds"))==F){

  cat("Classification binaire pour chaque classe \n") 
  folds_dt <- make_common_folds(
    ids = pam50$ID,
    k = 10,
    seed = 123
  )
  results <- pam50_rf_knn_classification(
  module_matrix_list = module_matrix_list,
  pam50 = pam50,
  folds_dt=folds_dt,
  k_knn = 20,
  classifier = "knn"
)
  setnames(results$formatted_table, old = "Method",new="Module Scoring Method")
  saveRDS(results, paste0(folder_for_res,"KNN_results.rds"))
  # write.csv2(results$formatted_table, paste0(folder_for_res, "knn_res_brca.csv"))
}else{
  results = readRDS(paste0(folder_for_res,"KNN_results.rds"))
  
}

# if(file.exists(paste0(folder_for_res,"multi_class_acc_KNN_results.rds"))==F){
#   
#   cat("Classification multiclasse \n")
#   folds_dt <- make_common_folds(
#     ids = pam50$ID,
#     k = 10,
#     seed = 123
#   )
# 
#   
#   accuracy_table <- pam50_rf_knn_accuracy_multiclass(
#     module_matrix_list = module_matrix_list,
#     pam50 = pam50,
#     folds_dt = folds_dt,
#     k_knn = 20,
#     classifier = "knn"
#   )
#   saveRDS(accuracy_table, paste0(folder_for_res,"multi_class_acc_KNN_results.rds"))
# }




# 9. Ordinality ------------------------------------------------------------

cat("compute ordinality analysis \n")

if(file.exists(paste0(folder_for_res, "ordinality_res.rds"))){
library(clinfun)
  
ordinality_results <- evaluate_pam50_ordinality_phate(
  module_matrix_list = module_matrix_list,
  pam50 = pam50
)
setnames(ordinality_results$results_table, old = "Method",new="Module Scoring Method")
saveRDS(ordinality_results, file=paste0(folder_for_res, "ordinality_res.rds"))
ggsave(plot = ordinality_results$combined_plot, filename=paste0(folder_for_res,"ordinality_comb.pdf"), height = 15, width = 15)
}else{
  ordinality_results = readRDS(paste0(folder_for_res, "ordinality_res.rds"))
}

#Fusion avec résultats KNN : 
results$formatted_table[ordinality_results$results_table, Kendall_tau :=i.Kendall_tau , on="Module Scoring Method" ]
results$formatted_table[ordinality_results$results_table, JT_statistic :=i.JT_statistic , on="Module Scoring Method" ]
results$formatted_table[ordinality_results$results_table, JT_pvalue :=i.JT_pvalue , on="Module Scoring Method" ]
results$formatted_table[,JT_pvalue := ifelse(
  JT_pvalue < 1e-16,
  "< 1e-16",
  formatC(JT_pvalue, format = "e", digits = 2)
)]
results$formatted_table[,Kendall_tau := round(Kendall_tau,digits = 3)]



write.csv2(results$formatted_table, paste0(folder_for_res, "knn_ordinality_res_brca.csv"))


cat("Finished ! \n")




# 10. Formatting figure ---------------------------------------------------

theme_cadre <- theme(
  panel.border = element_rect(
    colour = "black",
    fill   = NA,
    linewidth = 0.6
  ),
  plot.background = element_rect(
    colour = NA,
    fill   = NA
  ),
  panel.spacing = unit(2, "pt")
)



plots_umap_pam50 = make_umap_pam50(module_matrix_list = module_matrix_list,ids = ids, pam50 = pam50)








plots_ordinality = ordinality_results$plots_by_method
plots_ordinality$GEDO = plots_ordinality$GEDO+theme_cadre
plots_ordinality$GEDOcorr = plots_ordinality$GEDOcorr+theme_cadre
plots_ordinality$UMAP_GEDO = plots_ordinality$UMAP_GEDO+theme_cadre
plots_ordinality$PCA1 = plots_ordinality$PCA1+theme_cadre
plots_ordinality$MZS = plots_ordinality$MZS+theme_cadre
plots_ordinality$ssGSEA = plots_ordinality$ssGSEA+theme_cadre
plots_ordinality$GSVA = plots_ordinality$GSVA+theme_cadre





roc_curbes = results$roc_plots_by_class
roc_curbes$LumA = roc_curbes$LumA+theme_cadre
roc_curbes$LumB = roc_curbes$LumB+theme_cadre
roc_curbes$Her2 = roc_curbes$Her2+theme_cadre
roc_curbes$Basal = roc_curbes$Basal+theme_cadre


library(grid)
library(gridExtra)
library(ggplot2)

method_order=names(plots_ordinality)

grob_a <- arrangeGrob(
  grobs = plots_umap_pam50[method_order],
  ncol = 1,
  top = textGrob("a", x = 0, just = "left",
                 gp = gpar(fontface = "bold", fontsize = 14))
)

grob_b <- arrangeGrob(
  grobs = plots_ordinality[method_order],
  ncol = 1,
  top = textGrob("b", x = 0, just = "left",
                 gp = gpar(fontface = "bold", fontsize = 14))
)


grob_c <- arrangeGrob(
  grobs = roc_curbes[c("LumA","LumB","Her2","Basal")],
  ncol = 2,
  top = textGrob("c", x = 0, just = "left",
                 gp = gpar(fontface = "bold", fontsize = 14))
)

grob_ab <- arrangeGrob(
  grob_a,
  grob_b,
  ncol = 2
)



final_figure <- arrangeGrob(
  grob_ab,
  grob_c,
  ncol = 1,
  heights = c(3, 1)
)

# final_figure <- arrangeGrob(
#   grob_ab,
#   grob_c,
#   ncol = 2,widths = c(2,2)
# )



ggsave(
  filename = paste0(folder_for_res, "Figure_tcga_pam50.pdf"),
  plot = final_figure,
  width = 15,
  height = 30
)






