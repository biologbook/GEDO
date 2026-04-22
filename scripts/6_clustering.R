
# 1. Package and functions ------------------------------------------------


folder_script_to_source = "/home/clem/GEDO/R/"
source(file = paste0(folder_script_to_source,"GEDO.R"))
source(file = paste0(folder_script_to_source,"functions_article.R"))

detach_all_packages()
packages_list = c("data.table", "FNN", "magrittr", "tryCatchLog","rlist",
                  "igraph", "pbapply","rgl","rdist","msigdbr","msigdbdf",
                  "future.apply","dbscan", "progressr","pheatmap",
                  "gridExtra","MASS","plotly","ggplot2", "patchwork",
                  "caret","randomForest", "dplyr", "tidyr","parallel", "pROC", "ggpubr","fpc","cluster","phateR","reticulate", "devtools", "bench","RANN",
                  # "GSVA",
                  "BiocParallel", "visdat", "magick", "ggbreak","cowplot","grid","DBCVindex","clusterCrit")
sapply(packages_list, require, character.only = TRUE)


# 2. Data --------------------------------------------------



folder_for_res = "/home/clem/GEDO/results/"
num_cores=5




# 3. run gedo and others --------------------------------------------------
run_all_module_matrices = function(rna_seq_data, folder_for_res){
  
  
  message("Loading data")
  folder_to_data ="/home/clem/GEDO/data/"
  
  PS_brutes <- readRDS(paste0(folder_to_data,"PS_brutes.rds"))
  rna_seq_data = readRDS(paste0(folder_to_data, "bulk_rna_seq_final_data_batch_corrected_high_cv_scaled.rds"))
  
  rna_seq_data[, SAMPLING_OMIC_NUMBER:=paste0("N",SAMPLING_OMIC_NUMBER)]
  rna_seq_data[PS_brutes, diag := i.DIAGNOSIS_DISEASE_AT_ONSET, on="SAMPLING_OMIC_NUMBER"]
  rna_seq_data[PS_brutes, control:=i.DIAGNOSIS_ARM, on="SAMPLING_OMIC_NUMBER"]
  rna_seq_data[control=="Control", diag:="Control"]
  rna_seq_data=rna_seq_data[,diag:=factor(diag, levels = c("Control","SjS"))]
  rna_seq_data=rna_seq_data[diag %in% c("Control","SjS")]
  rna_seq_data[diag=="SjS",diag:="SjD"]
  diag = rna_seq_data$diag
  PS_brutes[,diag:=DIAGNOSIS_DISEASE_AT_ONSET]
  PS_brutes[DIAGNOSIS_ARM=="Control", diag:="Control"]
  diag=factor(diag,levels=c("SjD","Control"))
  omic_id = rna_seq_data$SAMPLING_OMIC_NUMBER
  ifn_score_table=data.table(SAMPLING_OMIC_NUMBER=omic_id)
  ifn_score_table[PS_brutes, IFN_score := i.EXPRESSION_PRECISESADS_IFN, on="SAMPLING_OMIC_NUMBER"]
  ifn_score=ifn_score_table$IFN_score
  rna_seq_data[, diag:=NULL][, SAMPLING_OMIC_NUMBER:=NULL][, control:=NULL]
  
  message("GEDO")
  
  
if(file.exists(paste0(folder_for_res, "gedo_obj.rds"))){
  gedo_obj = readRDS(paste0(folder_for_res, "gedo_obj.rds"))
}else{
  gedo_obj = gedo(data = rna_seq_data, diag = diag, reference_group = "Control",
                  category = "C7",subcategory = "IMMUNESIGDB",k_lof = 30,
                  core_pct = 0.9, k_graph = 15, dim_reduc_method="none", num_cores=num_cores,
                  distance="euclidean")
  
  saveRDS(gedo_obj, file=paste0(folder_for_res, "gedo_obj.rds"))
}
gedo_obj$diag=factor(diag,levels=c("SjD","Control"))


message("GEDOcorr")

if(file.exists(paste0(folder_for_res, "gedo_corr_obj.rds"))){
  gedo_corr_obj = readRDS(paste0(folder_for_res, "gedo_corr_obj.rds"))
}else{
  gedo_corr_obj = gedo(data = rna_seq_data, diag = diag, reference_group = "Control",
                       category = "C7",subcategory = "IMMUNESIGDB",k_lof = 30,
                       core_pct = 0.9, k_graph = 15, dim_reduc_method="none", num_cores=num_cores,
                       distance="correlation")
  
  saveRDS(gedo_corr_obj, file=paste0(folder_for_res, "gedo_corr_obj.rds"))
}
gedo_corr_obj$diag=factor(diag,levels=c("SjD","Control"))


message("UMAP_GEDO")


if(file.exists(paste0(folder_for_res, "umap_gedo_obj.rds"))){
  umap_gedo_obj = readRDS(paste0(folder_for_res, "umap_gedo_obj.rds"))
}else{
  umap_gedo_obj = gedo(data = rna_seq_data, diag = diag, reference_group = "Control",
                       category = "C7",subcategory = "IMMUNESIGDB",k_lof = 30,
                       core_pct = 0.9, k_graph = 15, dim_reduc_method="umap", ncomp = 10,
                       dim_reduc_dist_method = "correlation", num_cores=num_cores,
                       distance="euclidean")
  
  saveRDS(umap_gedo_obj, file=paste0(folder_for_res, "umap_gedo_obj.rds"))
}
umap_gedo_obj$diag=factor(diag,levels=c("SjD","Control"))


message("PCA1")

if(file.exists(paste0(folder_for_res, "pca1_obj.rds"))){
  pca1_module_matrix = readRDS(paste0(folder_for_res, "pca1_obj.rds"))
}else{
  pca1_module_matrix = compute_module_matrix(method="pca1", data=rna_seq_data, 
                                             diag=diag, 
                                             reference_group="Control",
                                             category="C7",
                                             subcategory = "IMMUNESIGDB",
                                             num_cores=num_cores)
  
  saveRDS(pca1_module_matrix, file=paste0(folder_for_res, "pca1_obj.rds"))
}


message("MZS")

if(file.exists(paste0(folder_for_res, "mzscore_obj.rds"))){
  mean_z_score_module_matrix=readRDS(paste0(folder_for_res, "mzscore_obj.rds"))
}else{
  mean_z_score_module_matrix = compute_module_matrix(method="mean_z_score", data=rna_seq_data, 
                                                     diag=diag, 
                                                     reference_group="Control",
                                                     category="C7",
                                                     subcategory = "IMMUNESIGDB",
                                                     num_cores=num_cores)
  
  saveRDS(mean_z_score_module_matrix, file=paste0(folder_for_res,"mzscore_obj.rds"))
}



message("ssGSEA")

expr_mat <- t(as.matrix(rna_seq_data))

immun_c7=data.table(msigdbr(species = "Homo sapiens", collection = "C7", subcollection = "IMMUNESIGDB"))
geneSets <- split(immun_c7$ensembl_gene, f = immun_c7$gs_name)

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
    config=list(param_ssgsea=param_ssgsea, reference_group="Control"),
    diag=as.factor(diag)
  )
  
  saveRDS(ssgsea_obj, file = paste0(folder_for_res, "ssgsea_obj.rds"))
}

message("GSVA")


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
  )
  
  scores_gsva  <- gsva(param_gsva, BPPARAM=MulticoreParam(workers = num_cores,progressbar = T))
  
  scores_gsva=data.table(t(scores_gsva))
  scores_gsva[, (names(scores_gsva)) := lapply(.SD, function(x) (x - min(x, na.rm=TRUE)) / (max(x, na.rm=TRUE) - min(x, na.rm=TRUE)))]
  
  gsva_obj=list(
    module_matrix = scores_gsva,
    config=list(param_gsva=param_gsva, reference_group="Control"),
    diag=as.factor(diag)
  )
  
  saveRDS(gsva_obj, file = paste0(folder_for_res, "gsva_obj.rds"))
}

rm(rna_seq_data)

message("Saving module matrices")

module_matrix_list = list(GEDO=gedo_obj,
                 
                   GEDOcorr = gedo_corr_obj,
                   
                   UMAP_GEDO=umap_gedo_obj,
                   
                   PCA1=pca1_module_matrix,
                   
                   MEAN_Z_SCORES=mean_z_score_module_matrix,
                   
                   ssGSEA = ssgsea_obj,
                   
                   GSVA = gsva_obj
                   )

saveRDS(module_matrix_list, file = paste0(folder_for_res, "module_matrix_list.rds"))

}



if(file.exists(paste0(folder_for_res, "module_matrix_list.rds"))){
  module_matrix_list = readRDS(paste0(folder_for_res, "module_matrix_list.rds"))
}else{
run_all_module_matrices(rna_seq_data=rna_seq_data, folder_for_res=folder_for_res)
cat("Finished !")
module_matrix_list = readRDS(paste0(folder_for_res, "module_matrix_list.rds"))
}


# 4. PCA and clustering ---------------------------------------------------


  
#' Run clustering and plot clusters in PCA embeddings
#' @param pca_df PCA embedding
#' @param mat matrix to clusterise
#' @param k number of clusters
#' @param element_name names of methods to test
plot_pca_cluster_2 <- function(pca_df,mat, k, element_name) {
  
  
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
  
  # Clustering hiérarchique
  message("c=", k, "step : Clustering")
  hc <- hclust(dist(mat))
  pca_df$cluster <- factor(cutree(hc, k = k))
  
  # Plot
  ggplot(pca_df, aes(PC1, PC2, color = cluster)) +
    geom_point(size = 2, alpha = 0.5) +
    labs(
      title = paste0("c = ", k),
      x = "PC1",
      y = "PC2",
      color = "Cluster"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 10),
      legend.position = "none"
    ) + theme_cadre + 
    theme(legend.position = "none")
}


stopifnot(length(module_matrix_list) == 7)
stopifnot(!is.null(names(module_matrix_list)))

all_rows <- list()

for (element_name in names(module_matrix_list)) {
  print(element_name)
  mat <- module_matrix_list[[element_name]]$module_matrix


  mat <- as.matrix(mat)

  #PCA
  message("PCA")
  set.seed(123)
  pca <- prcomp(mat, scale. = TRUE)
  pca_df <- data.frame(
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2]
  )

  # k = 2 to 5
  plots_k <- lapply(2:5, function(k) {
    plot_pca_cluster_2(pca_df, mat, k, element_name)
  })

 
  row_plot <- wrap_plots(plots_k, ncol = 4)
  all_rows[[element_name]] <- row_plot
}

# 
# final_plot <- wrap_plots(all_rows, ncol = 1)
# 
# final_plot
# 
# ggsave(plot=final_plot,filename=paste0(folder_for_res, "pca_clusters.pdf"),width = 30, height=30 )
# # }



#' Prepare a row for patchwork
#' @param plot_list list of plots
#' @param label label for the row
make_row <- function(plot_list, label) {
  
  label_grob <- grid::textGrob(
    label,
    x = 0, just = "left",
    gp = grid::gpar(fontsize = 14, fontface = "bold")
  )
  
  wrap_plots(
    wrap_elements(label_grob),
    plot_list[[1]],
    plot_list[[2]],
    plot_list[[3]],
    plot_list[[4]],
    
    ncol = 5,
    widths = c(0.15, 0.44, 0.44, 0.44, 0.44)
  )
}




names(all_rows) = gsub(x= names(all_rows), pattern="MEAN_Z_SCORES",replacement="MZS")
rows <- Map(
  make_row,
  plot_list = all_rows,
  label     = names(all_rows)
)

final_plot <- wrap_plots(rows, ncol = 1) +
  plot_layout(guides = "collect") &
  theme(legend.position = "none")

ggsave(plot=final_plot, filename = paste0(folder_for_res,"pca_clusters.pdf"), width=20, height = 20)








# 5. Clustering evaluation ------------------------------------------------

library(DBCVindex)


if(file.exists(paste0(folder_for_res,"clustering_quality.rds"))){
  cl_res = readRDS(paste0(folder_for_res,"clustering_quality.rds"))
}else{
  cl_res = compute_clustering_quality(matrix_list = module_matrix_list, num_cores=num_cores)
  saveRDS(cl_res, file=paste0(folder_for_res,"clustering_quality.rds"))
}



# 6. Saving plot ----------------------------------------------------------


ggsave(plot=cl_res$plot, filename=paste0(folder_for_res, "clustering_quality_res.pdf"), width = 15, height=8)



