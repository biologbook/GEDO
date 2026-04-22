
# 1. Packages and scripts -------------------------------------------------

folder_script_to_source = "/home/clem/GEDO/R/GEDO/"
source(file = paste0(folder_script_to_source,"GEDO.R"))
source(file = paste0(folder_script_to_source,"functions_article.R"))

detach_all_packages()
packages_list = c("data.table", "FNN", "magrittr", "tryCatchLog","rlist",
                  "igraph", "pbapply","rgl","rdist","msigdbr","msigdbdf",
                  "future.apply","dbscan", "progressr","pheatmap",
                  "gridExtra","MASS","plotly","ggplot2", "patchwork",
                  "caret","randomForest", "dplyr", "tidyr","parallel",
                  "pROC", "ggpubr","fpc","cluster","phateR","reticulate",
                  "devtools", "bench","RANN", "GSVA", "BiocParallel",
                  "visdat", "magick", "ggbreak","cowplot","grid","progress")

install_and_load(packages_list)


# 2. Data -----------------------------------------------------------------

folder_for_res = "/home/clem/GEDO/article_GEDO/results/"
folder_to_data ="/home/clem/GEDO/article_GEDO/data/"

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


num_cores=3





# 3. Functions ------------------------------------------------------------


generate_gedo_configs <- function() {
  expand.grid(
    method     = c("GEDO", "GEDOcorr"),
    core_pct   = c(1.0, 0.9, 0.8, 0.7),
    k_lof      = c(5, 10, 20, 30),
    k_graph    = c(5, 10, 20, 30),
    stringsAsFactors = FALSE
  )
}

#' Run GEDO on parameter grid
#' @param rna_seq_data clinical data to test
#' @param diag vector with diagnosis (SjD or Control)
#' @param output_dir directory for output
#' @param category GSEA category
#' @param subcategory GSEA subcategory
#' @param dim_reduc_method dimension reduction method : "none" or "umap"
#' @param num_cores number of cores for parallelization
#' @param overwrite overwrite files when saving
#' @param umap TRUE/FALSE. If TRUE, run UMAP

run_gedo_grid <- function(
    rna_seq_data,
    diag,
    output_dir,
    category = "C7",
    subcategory = "IMMUNESIGDB",
    dim_reduc_method = "none",
    num_cores = 1,
    overwrite = FALSE,
    umap=F
) {
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  if(umap==F){
  configs <- generate_gedo_configs()
  }
  
  if(umap==T){
    configs=generate_umap_gedo_configs()
  }
  
  n_runs <- nrow(configs)
  
  pb <- progress::progress_bar$new(
    format = "  [:bar] :percent | :current/:total | elapsed: :elapsed | eta: :eta",
    total = n_runs,
    clear = FALSE,
    width = 60
  )
  
  results <- vector("list", nrow(configs))
  
  for (i in seq_len(nrow(configs))) {
    
    cfg <- configs[i, ]
    pb$tick(tokens = list(
      cfg = paste(cfg$method,
                  "core", cfg$core_pct,
                  "klof", cfg$k_lof,
                  "kknn", cfg$k_graph)
    ))
    
    
    
    
    distance <- ifelse(cfg$method == "GEDO",
                       "euclidean",
                       "correlation")
    
    file_name <- sprintf(
      "%s_core%s_klof%s_kknn%s.rds",
      cfg$method,
      cfg$core_pct,
      cfg$k_lof,
      cfg$k_graph
    )
    
    file_path <- file.path(output_dir, file_name)
    
    if (file.exists(file_path) && !overwrite) {
      message("[SKIP] ", file_name)
      next
    }
    
    message("[RUN] ", file_name)
    
    gedo_obj <- gedo(
      data = rna_seq_data,
      diag = diag,
      reference_group = "Control",
      category = category,
      subcategory = subcategory,
      k_lof = cfg$k_lof,
      core_pct = cfg$core_pct,
      k_graph = cfg$k_graph,
      dim_reduc_method = dim_reduc_method,
      distance = distance,
      num_cores = num_cores
    )
    
    saveRDS(gedo_obj, file = file_path)
    
    results[[i]] <- list(
      config = cfg,
      file = file_path
    )
  }
  
  invisible(results)
}




# 4. Run function for GEDO and GEDOcorr---------------------------------------------------------


set.seed(1234)


run_gedo_grid(
  rna_seq_data = rna_seq_data,
  diag = diag,
  output_dir = folder_for_res,
  num_cores = num_cores
)



# 5. Run functions for UMAP_GEDO ------------------------------------------

generate_umap_gedo_configs <- function() {
  expand.grid(
    method     = "UMAP_GEDO",
    core_pct   = 0.9,
    k_lof      = 30,
    k_graph    = c(5, 10, 20, 30),
    umap_ncomp = c(2,10,20,30),
    stringsAsFactors = FALSE
  )
}

#' Run UMAP-GEDO on parameter grid
#' @param rna_seq_data clinical data to test
#' @param diag vector with diagnosis (SjD or Control)
#' @param output_dir directory for output
#' @param category GSEA category
#' @param subcategory GSEA subcategory
#' @param dim_reduc_method dimension reduction method : "none" or "umap"
#' @param num_cores number of cores for parallelization
#' @param overwrite overwrite files when saving
#' @param umap TRUE/FALSE. If TRUE, run UMAP
run_umap_gedo_grid <- function(
    rna_seq_data,
    diag,
    output_dir,
    category = "C7",
    subcategory = "IMMUNESIGDB",
    dim_reduc_method = "none",
    num_cores = 1,
    overwrite = FALSE,
    umap = FALSE
) {
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  if (!umap) {
    configs <- generate_gedo_configs()
  } else {
    configs <- generate_umap_gedo_configs()
    dim_reduc_method <- "umap"
  }
  
  n_runs <- nrow(configs)
  
  pb <- progress::progress_bar$new(
    format = "[:bar] :percent | :current/:total | elapsed: :elapsed | eta: :eta",
    total  = n_runs,
    clear  = FALSE,
    width  = 60
  )
  
  results <- vector("list", n_runs)
  
  for (i in seq_len(n_runs)) {
    
    cfg <- configs[i, ]
    
    pb$tick(tokens = list(
      cfg = paste0(
        cfg$method,
        "_kknn", cfg$k_graph,
        if (umap) paste0("_ncomp", cfg$ncomp) else ""
      )
    ))
    
    distance <- "euclidean"
    
    ## ---- file naming (cohérent avec l’existant)
    if (!umap) {
      file_name <- sprintf(
        "%s_core%s_klof%s_kknn%s.rds",
        cfg$method,
        cfg$core_pct,
        cfg$k_lof,
        cfg$k_graph
      )
    } else {
      file_name <- sprintf(
        "%s_core%s_klof%s_kknn%s_ncomp%s.rds",
        cfg$method,
        cfg$core_pct,
        cfg$k_lof,
        cfg$k_graph,
        cfg$umap_ncomp
      )
    }
    
    file_path <- file.path(output_dir, file_name)
    
    if (file.exists(file_path) && !overwrite) {
      message("[SKIP] ", file_name)
      next
    }
    
    message("[RUN] ", file_name)
    
    gedo_obj <- gedo(
      data = rna_seq_data,
      diag = diag,
      reference_group = "Control",
      category = category,
      subcategory = subcategory,
      distance = distance,
      dim_reduc_method = dim_reduc_method,
      dim_reduc_dist_method="correlation",
      ncomp = if (umap) cfg$umap_ncomp else NULL,
      k_lof = cfg$k_lof,
      core_pct = cfg$core_pct,
      k_graph = cfg$k_graph,
      num_cores = num_cores
    )
    
    saveRDS(gedo_obj, file = file_path)
    
    results[[i]] <- list(
      config = cfg,
      file   = file_path
    )
  }
  
  invisible(results)
}

            




run_umap_gedo_grid(
  rna_seq_data = rna_seq_data,
  diag = diag,
  output_dir = folder_for_res,
  num_cores = num_cores,
  umap=T,
  overwrite = FALSE
)
