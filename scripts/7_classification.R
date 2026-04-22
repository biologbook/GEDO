
# 1. packages and data ----------------------------------------------------

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
                  "BiocParallel", "visdat", "magick", "ggbreak","cowplot","grid", "DBCVindex","clusterCrit")

print(sapply(packages_list, require, character.only = TRUE))

matrix_list_raw = readRDS("/home/clem/GEDO/results/matrix_list_without_noize.rds")


if(file.exists("/home/clem/GEDO/results/meta_data.rds")){
  meta_data = readRDS("/home/clem/GEDO/results/meta_data.rds")
}else{

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
# saveRDS(diag, "/shared/projects/toposads/finalresult/7_gedo_reviewing/figures/diag.rds")
PS_brutes[,diag:=DIAGNOSIS_DISEASE_AT_ONSET]
PS_brutes[DIAGNOSIS_ARM=="Control", diag:="Control"]
diag=factor(diag,levels=c("SjD","Control"))
omic_id = rna_seq_data$SAMPLING_OMIC_NUMBER
# saveRDS(omic_id, "/shared/projects/toposads/finalresult/7_gedo_reviewing/figures/omic_id.rds")

meta_data=data.table(SAMPLING_OMIC_NUMBER=omic_id)
meta_data$diag = diag
meta_data[PS_brutes, IFN_score := i.EXPRESSION_PRECISESADS_IFN, on="SAMPLING_OMIC_NUMBER"]
meta_data[PS_brutes, SSA_52 := i.AUTOANTIBODY_SSA_52_CALL, on="SAMPLING_OMIC_NUMBER"]
meta_data[PS_brutes, SSA_60 := i.AUTOANTIBODY_SSA_60_CALL, on="SAMPLING_OMIC_NUMBER"]
meta_data[PS_brutes, SSA := i.AUTOANTIBODY_SSA_CALL, on="SAMPLING_OMIC_NUMBER"]
meta_data[PS_brutes, SSB := i.AUTOANTIBODY_SSB_CALL, on="SAMPLING_OMIC_NUMBER"]
meta_data[PS_brutes, ESSDAI := i.SJS_ESSDAI, on="SAMPLING_OMIC_NUMBER"]

saveRDS(meta_data, "/home/clem/GEDO/results/meta_data.rds")
}

# folder_for_res = "/shared/projects/toposads/finalresult/7_gedo_reviewing/classification/"
folder_for_res = "/home/clem/GEDO/results/"
num_cores=5




# 2. Functions ------------------------------------------------------------


#' Run KNN and RF classification
#' @param model "knn" or "rf"
#' @param k number of neighbors for KNN
#' @param k_folds number of folds
#' @param dt_list list of module matrices
#' @param num_cores number of cores for parallelization
#' @param y class to predict
compute_classification_with_ci <- function(model = "knn", k, k_folds = 10, dt_list, num_cores = NULL, y) {
  # Set default number of cores if not specified
  if (is.null(num_cores)) {
    num_cores <- detectCores() - 1
  }
  
  print("computing k-folds")
  
  roc_curves <- lapply(names(dt_list), function(dt_name) {
    cat(paste0(dt_name, "\n"))
    
    # Load the dataset and add target to the feature matrix
    
    dt <- data.table::copy(dt_list[[dt_name]]$module_matrix)
    stopifnot(nrow(dt) == length(y))

    
    # Attach target as column 'y'
    if ("y" %in% names(dt)) {
      message("y column was already there")
      dt[, y := NULL]
    }
    
    dt[, y :=as.factor(y)]
    dt = dt[!is.na(y)]
    
    # Split 70% train / 30% test
    set.seed(423)
    trainIndex <- createDataPartition(dt$y, p = 0.7, list = FALSE)
    trainData <- dt[trainIndex, ]
    testData <- dt[-trainIndex, ]
    
    # Create k folds from the training data
    set.seed(423)
    folds <- createFolds(trainData$y, k = k_folds, list = TRUE)
    
    all_probs <- c()
    all_labels <- c()
    
    # Setup parallel processing
    plan(multisession, workers = num_cores)
    #or
    # plan(sequential)
    
    handlers("txtprogressbar")
    
    cat(paste0("Metric : ", dt_name), "\n")
    cat("Parallelization on folds", "\n")
    
    export = list(model = model, k = k, k_folds = k_folds, folds = folds, trainData = trainData)
    
    # Run training/testing in parallel across folds
    fold_results <- with_progress({
      p <- progressor(along = folds)
      future_lapply(1:length(folds), function(i) {
        cat(paste0("Fold : ", i), "\n")
        
        testIndexes <- folds[[i]]
        testFold <- trainData[testIndexes, ]
        trainFold <- trainData[-testIndexes, ]
        
        # Train model
        if (model == "knn") {
          classifier <- knn3(y ~ ., data = trainFold, k = k)
        } else if (model == "rf") {
          classifier <- randomForest(as.formula(y ~ .), data = trainFold, ntree = 400, probability = TRUE)
        } else {
          stop("Unsupported model type. Choose 'knn' or 'rf'.")
        }
        
        # Predict probabilities on test fold
        probs <- predict(classifier, newdata = testFold, type = "prob")[, 2]
        return(list(probs = probs, labels = testFold$y))
      }, future.seed = TRUE, future.globals = export, future.packages = c("randomForest", "caret"))
    })
    
    # Aggregate probabilities and labels across folds
    all_probs <- unlist(lapply(fold_results, `[[`, "probs"))
    all_labels <- unlist(lapply(fold_results, `[[`, "labels"))
    
    # Compute ROC and AUC
    roc_curve <- roc(all_labels, all_probs)
    
    # Compute 95% Confidence Interval of AUC
    auc_ci <- ci.auc(roc_curve)
    cat(paste0("AUC for ", dt_name, " : ", round(auc(roc_curve), 3),
               " [", round(auc_ci[1], 3), " - ", round(auc_ci[3], 3), "]\n"))
    
    # Attach CI to ROC object
    roc_curve$ci <- auc_ci
    return(roc_curve)
  })
  
  names(roc_curves) <- names(dt_list)
  print("roc_curves_ok")
  
  # Summary table of AUC and CI for each dataset
  res_table <- data.table(
    methods = names(dt_list),
    AUC = as.numeric(NA),
    CI_lower = as.numeric(NA),
    CI_upper = as.numeric(NA)
  )
  for (method in names(roc_curves)) {
    auc <- as.numeric(roc_curves[[method]]$auc)
    ci_vals <- as.vector(roc_curves[[method]]$ci)
    res_table[methods == method, `:=`(AUC = auc, CI_lower = ci_vals[1], CI_upper = ci_vals[3])]
  }
  setorder(res_table, -AUC)
  
  # Prepare ROC curves for plotting
  print("printing roc curves")
  roc_data <- do.call(rbind, mclapply(names(roc_curves), function(dt_name) {
    data.frame(
      inv_specificity = 1 - roc_curves[[dt_name]]$specificities,
      sensibility = roc_curves[[dt_name]]$sensitivities,
      Dataset = dt_name
    )
  }, mc.cores = num_cores))
  
  roc_data <- data.table(roc_data)
  desired_order <- names(dt_list)
  roc_data[, Dataset := factor(Dataset, levels = desired_order)]
  
  # Plot ROC curves
  if (model == "rf") {
    plot <- ggplot(roc_data, aes(x = inv_specificity, y = sensibility, color = Dataset)) +
      geom_line(size = 0.5) +
      scale_colour_viridis_d(option = "viridis")+
      labs(x = "1 - Specificity", y = "Sensitivity", color = "Module Scoring Method") +
      theme_minimal() +
      guides(color = "none")
  }
  if (model == "knn") {
    plot <- ggplot(roc_data, aes(x = inv_specificity, y = sensibility, color = Dataset)) +
      geom_line(size = 0.5) +
      scale_colour_viridis_d(option = "viridis")+
      labs(x = "1 - Specificity", y = "Sensitivity", color = "Module Scoring Method") +
      theme_minimal()
    
  }
  
  # Return summary results and plot
  list_to_return = list(auc_results = res_table, roc_curves = plot, roc_curves_obj=roc_curves)
  return(list_to_return)
}

#' Run KNN and RF regression of ifn score
#' @param model "knn" or "rf"
#' @param k number of neighbors for KNN
#' @param k_folds number of folds
#' @param dt_list list of module matrices
#' @param num_cores number of cores for parallelization
#' @param y class to predict
compute_regression_with_cv <- function(
    model = "knn",
    k,
    k_folds = 50,
    dt_list,
    num_cores = NULL,
    y
) {
  
  library(data.table)
  library(caret)
  library(randomForest)
  
  if (is.null(num_cores)) {
    num_cores <- 1
  }
  
  message("Computing k-fold CV (regression)")
  
  res_list <- lapply(names(dt_list), function(dt_name) {
    
    cat(paste0("Method: ", dt_name, "\n"))
    
    # ---- Load & protect data ----
    dt <- data.table::copy(dt_list[[dt_name]]$module_matrix)
    
    # Sanity checks
    stopifnot(nrow(dt) == length(y))
    
    # Ensure no leftover target
    if ("y" %in% names(dt)) {
      dt[, y := NULL]
    }
    
    # Attach target
    dt[, y := as.numeric(y)]
    dt <- dt[!is.na(y)]
    
    # ---- Train / test split (train only used for CV) ----
    set.seed(423)
    train_idx <- caret::createDataPartition(dt$y, p = 0.7, list = FALSE)
    trainData <- dt[train_idx, ]
    
    # ---- K-fold CV on training set ----
    set.seed(423)
    folds <- caret::createFolds(trainData$y, k = k_folds, list = TRUE)
    
    fold_results <- lapply(seq_along(folds), function(i) {
      
      test_idx  <- folds[[i]]
      testFold <- trainData[test_idx, ]
      trainFold <- trainData[-test_idx, ]
      
      # ---- Train model ----
      if (model == "knn") {
        
        regressor <- caret::knnreg(
          x = trainFold[, !"y", with = FALSE],
          y = trainFold$y,
          k = k
        )
        
        preds <- predict(
          regressor,
          newdata = testFold[, !"y", with = FALSE]
        )
        
      } else if (model == "rf") {
        
        regressor <- randomForest::randomForest(
          y ~ .,
          data = trainFold,
          ntree = 400
        )
        
        preds <- predict(regressor, newdata = testFold)
        
      } else {
        stop("Unsupported model type. Choose 'knn' or 'rf'.")
      }
      
      obs <- testFold$y
      
      # ---- Metrics ----
      rmse <- sqrt(mean((obs - preds)^2))
      mae  <- mean(abs(obs - preds))
      
      list(
        fold = i,
        rmse = rmse,
        mae  = mae
      )
    })
    
    # ---- Fold-level table ----
    folds_dt <- rbindlist(lapply(fold_results, function(res) {
      data.table(
        method = dt_name,
        fold   = res$fold,
        RMSE   = res$rmse,
        MAE    = res$mae
      )
    }))
    
    # ---- Summary table ----
    summary_dt <- folds_dt[, .(
      RMSE_mean = mean(RMSE),
      RMSE_sd   = sd(RMSE),
      MAE_mean  = mean(MAE),
      MAE_sd    = sd(MAE)
    ), by = method]
    
    list(
      summary = summary_dt,
      folds   = folds_dt
    )
  })
  
  # ---- Aggregate across methods ----
  summary_table <- rbindlist(lapply(res_list, `[[`, "summary"))
  folds_table   <- rbindlist(lapply(res_list, `[[`, "folds"))
  
  setorder(summary_table, RMSE_mean)
  
  message("Regression CV completed")
  
  return(list(
    summary = summary_table,
    folds   = folds_table
  ))
}


    


# 2. UMAP -----------------------------------------------------------------


cat("Computing umaps \n")

if(!file.exists(paste0(folder_for_res, "matrix_list.rds"))){
matrix_list<- lapply(matrix_list_raw, function(x) {
  set.seed(123)
  umap = data.table(uwot::umap(X = x$module_matrix,n_components = 10,n_neighbors = 10, metric = "euclidean"))
  x$module_matrix = umap
  return(x)
})
saveRDS(matrix_list, file = paste0(folder_for_res, "matrix_list.rds"))
}else{
  matrix_list=readRDS(paste0(folder_for_res, "matrix_list.rds"))
}


if(!file.exists(paste0(folder_for_res, "matrix_list_sjd.rds"))){
  matrix_list_sjd <- lapply(matrix_list_raw, function(x) {
    idx <- x$diag == "SjD"
    module_matrix <- x$module_matrix[idx, ]
    set.seed(123)
    umap = data.table(uwot::umap(X = module_matrix,n_components = 10,n_neighbors = 10, metric = "euclidean"))
    x$module_matrix = umap
    x$diag <- x$diag[idx]
    return(x)
  })
}else{
  matrix_list_sjd=readRDS(paste0(folder_for_res, "matrix_list_sjd.rds"))
}





# 3. Categorical Y --------------------------------------------------------



# 3a. KNN CTRL + SjD-----------------------------------------------------------------
k_ctrl_sjd = round(sqrt(length(meta_data$diag)))

if(!file.exists(paste0(folder_for_res, "knn_ssa52.rds"))){
knn_ssa52 = compute_classification_with_ci(model = "knn", k=k_ctrl_sjd, k_folds = 10, dt_list=matrix_list, num_cores=num_cores, y=meta_data$SSA_52)
saveRDS(knn_ssa52, file=paste0(folder_for_res, "knn_ssa52.rds"))
}else{
  knn_ssa52=readRDS(paste0(folder_for_res, "knn_ssa52.rds"))
}

if(!file.exists(paste0(folder_for_res, "knn_ssa60.rds"))){
  knn_ssa60 = compute_classification_with_ci(model = "knn", k=k_ctrl_sjd, k_folds = 10, dt_list=matrix_list, num_cores=num_cores, y=meta_data$SSA_60)
  saveRDS(knn_ssa60, file=paste0(folder_for_res, "knn_ssa60.rds"))
}else{
  knn_ssa60=readRDS(paste0(folder_for_res, "knn_ssa60.rds"))
} 


if(!file.exists(paste0(folder_for_res, "knn_ssa.rds"))){
  knn_ssa = compute_classification_with_ci(model = "knn", k=k_ctrl_sjd, k_folds = 10, dt_list=matrix_list, num_cores=num_cores, y=meta_data$SSA)
  saveRDS(knn_ssa, file=paste0(folder_for_res, "knn_ssa.rds"))
}else{
  knn_ssa=readRDS(paste0(folder_for_res, "knn_ssa.rds"))
}

if(!file.exists(paste0(folder_for_res, "knn_ssb.rds"))){
  knn_ssb = compute_classification_with_ci(model = "knn", k=k_ctrl_sjd, k_folds = 10, dt_list=matrix_list, num_cores=num_cores, y=meta_data$SSB)
  saveRDS(knn_ssb, file=paste0(folder_for_res, "knn_ssb.rds"))
}else{
  knn_ssb=readRDS(paste0(folder_for_res, "knn_ssb.rds"))
}


if(!file.exists(paste0(folder_for_res, "knn_ifn.rds"))){
  knn_ifn = compute_regression_with_cv(model = "knn", k=k_ctrl_sjd, k_folds = 50, dt_list=matrix_list, num_cores=num_cores, y=meta_data$IFN_score)
  saveRDS(knn_ifn, file=paste0(folder_for_res, "knn_ifn.rds"))
}else{
  knn_ifn=readRDS(paste0(folder_for_res, "knn_ifn.rds"))
}

if(!file.exists(paste0(folder_for_res, "knn_essdai.rds"))){
  knn_essdai = compute_regression_with_cv(model = "knn", k=k_ctrl_sjd, k_folds = 50, dt_list=matrix_list, num_cores=num_cores, y=meta_data$ESSDAI)
  saveRDS(knn_essdai, file=paste0(folder_for_res, "knn_essdai.rds"))
}else{
  knn_essdai=readRDS(paste0(folder_for_res, "knn_essdai.rds"))
}




# 3b. KNN SjD only --------------------------------------------------------

k_sjd = round(sqrt(nrow(meta_data[diag=="SjD"])))

idx_sjd <- meta_data$diag == "SjD"


# matrix_list_sjd <- lapply(matrix_list, function(x) {
#   idx <- x$diag == "SjD"
#   x$module_matrix <- x$module_matrix[idx, ]
#   x$diag <- x$diag[idx]
#   x
# })
# matrix_list_sjd <- lapply(matrix_list, function(x) {
#   x$module_matrix <- x$module_matrix[idx_sjd, ]
#   x$diag <- x$diag[idx_sjd]
#   x
# })

y_ifn_sjd    <- meta_data$IFN_score[idx_sjd]
y_essdai_sjd <- meta_data$ESSDAI[idx_sjd]



if(!file.exists(paste0(folder_for_res, "knn_ssa52_sjd.rds"))){
  knn_ssa52_sjd = compute_classification_with_ci(model = "knn", k=k_sjd, k_folds = 10, dt_list=matrix_list_sjd, num_cores=num_cores, y=meta_data[diag=="SjD"]$SSA_52)
  saveRDS(knn_ssa52_sjd, file=paste0(folder_for_res, "knn_ssa52_sjd.rds"))
}else{
  knn_ssa52_sjd=readRDS(paste0(folder_for_res, "knn_ssa52_sjd.rds"))
}

if(!file.exists(paste0(folder_for_res, "knn_ssa60_sjd.rds"))){
  knn_ssa60_sjd = compute_classification_with_ci(model = "knn", k=k_sjd, k_folds = 10, dt_list=matrix_list_sjd, num_cores=num_cores, y=meta_data[diag=="SjD"]$SSA_60)
  saveRDS(knn_ssa60_sjd, file=paste0(folder_for_res, "knn_ssa60_sjd.rds"))
}else{
  knn_ssa60_sjd=readRDS(paste0(folder_for_res, "knn_ssa60_sjd.rds"))
} 


if(!file.exists(paste0(folder_for_res, "knn_ssa_sjd.rds"))){
  knn_ssa_sjd = compute_classification_with_ci(model = "knn", k=k_sjd, k_folds = 10, dt_list=matrix_list_sjd, num_cores=num_cores, y=meta_data[diag=="SjD"]$SSA)
  saveRDS(knn_ssa_sjd, file=paste0(folder_for_res, "knn_ssa_sjd.rds"))
}else{
  knn_ssa_sjd=readRDS(paste0(folder_for_res, "knn_ssa_sjd.rds"))
}

if(!file.exists(paste0(folder_for_res, "knn_ssb_sjd.rds"))){
  knn_ssb_sjd = compute_classification_with_ci(model = "knn", k=k_sjd, k_folds = 10, dt_list=matrix_list_sjd, num_cores=num_cores, y=meta_data[diag=="SjD"]$SSB)
  saveRDS(knn_ssb_sjd, file=paste0(folder_for_res, "knn_ssb_sjd.rds"))
}else{
  knn_ssb_sjd=readRDS(paste0(folder_for_res, "knn_ssb_sjd.rds"))
}

if(!file.exists(paste0(folder_for_res, "knn_ifn_sjd.rds"))){
  stopifnot(
    nrow(matrix_list_sjd[[1]]$module_matrix) ==
      sum(meta_data$diag == "SjD")
  )
  
  knn_ifn_sjd = compute_regression_with_cv(model = "knn", k=k_sjd, k_folds = 50, dt_list=matrix_list_sjd, num_cores=num_cores, y= y_ifn_sjd)
  saveRDS(knn_ifn_sjd, file=paste0(folder_for_res, "knn_ifn_sjd.rds"))
}else{
  knn_ifn_sjd=readRDS(paste0(folder_for_res, "knn_ifn_sjd.rds"))
}


if(!file.exists(paste0(folder_for_res, "knn_essdai_sjd.rds"))){
  stopifnot(
    nrow(matrix_list_sjd[[1]]$module_matrix) ==
      sum(meta_data$diag == "SjD")
  )
  knn_essdai_sjd = compute_regression_with_cv(model = "knn", k=k_sjd, k_folds = 50, dt_list=matrix_list_sjd, num_cores=num_cores, y=y_essdai_sjd)
  saveRDS(knn_essdai_sjd, file=paste0(folder_for_res, "knn_essdai_sjd.rds"))
}else{
  knn_essdai_sjd=readRDS(paste0(folder_for_res, "knn_essdai_sjd.rds"))
}


cat("KNN Finished ! \n")


# 3c KNN. Figure ---------------------------------------------------------------


a <- wrap_elements(full = as_grob(knn_ssa52$roc_curves))
b <- wrap_elements(full = as_grob(knn_ssa60$roc_curves))
c <- wrap_elements(full = as_grob(knn_ssa$roc_curves))
d <- wrap_elements(full = as_grob(knn_ssb$roc_curves))
e <- wrap_elements(full = as_grob(knn_ssa52_sjd$roc_curves))
f <- wrap_elements(full = as_grob(knn_ssa60_sjd$roc_curves))
g <- wrap_elements(full = as_grob(knn_ssa_sjd$roc_curves))
h <- wrap_elements(full = as_grob(knn_ssb_sjd$roc_curves))

comb <- a / b / c / d /e /f / g/ h     +
  plot_layout(ncol = 4, guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(
    plot.tag          = element_text(size = 20, face = "bold"),
    plot.tag.position = c(0, 1)
  )
ggsave(plot = comb, filename=paste0(folder_for_res, "knn_roc_curves_ab.pdf"), height=10, width = 25)




# # 4. RF -------------------------------------------------------------------
# 
# # 4a. RF CTRL + SjD-----------------------------------------------------------------

#--------
ordre_method <- knn_ifn$folds %>%
  group_by(method) %>%
  summarise(med_RMSE = median(RMSE, na.rm = TRUE)) %>%
  arrange(med_RMSE) %>%
  pull(method)

knn_ifn$folds$method <- factor(knn_ifn$folds$method, levels = ordre_method)

ordre_couleur <- unique(knn_ifn$folds$method)

rmse_ifn=ggplot(knn_ifn$folds, aes(x = method, y = RMSE, color = method)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.6) +
  scale_colour_viridis_d(
    option = "viridis",
    limits = ordre_couleur
  ) +
  theme_minimal()



rmse_ifn_sjd=ggplot(knn_ifn_sjd$folds, aes(x = method, y = RMSE, color = method)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.6) +
  scale_colour_viridis_d(
    option = "viridis",
    limits = ordre_couleur
  ) +
  theme_minimal()


library(dplyr)
library(ggplot2)
library(rlang)


#' Plot KNN performance metrics
#' @param data knn results
#' @param metric name of the metric
#' @param dataset_name name of the method tested
plot_knn_metric <- function(data, metric, dataset_name) {
  
  metric_sym <- sym(metric)
  
  df <- data$folds
  
  # ordre des méthodes par médiane de la métrique
  ordre_method <- df %>%
    group_by(method) %>%
    summarise(med = median(!!metric_sym, na.rm = TRUE)) %>%
    arrange(med) %>%
    pull(method)
  
  # ordre couleur original
  ordre_couleur <- unique(df$method)
  
  df <- df %>%
    mutate(
      method = factor(method, levels = ordre_method),
      metric_value = !!metric_sym,
      metric = metric,
      dataset = dataset_name
    )
  
  ggplot(df, aes(method, metric_value, color = method)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.6) +
    scale_colour_viridis_d(
      option = "viridis",
      limits = ordre_couleur
    ) +
    labs(
      x = NULL,
      y = metric,
      title = paste(dataset_name, "-", metric)
    ) +
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
}


p1 <- plot_knn_metric(knn_ifn,     "RMSE", "IFN")
p2 <- plot_knn_metric(knn_ifn,     "MAPE", "IFN")
p3 <- plot_knn_metric(knn_ifn_sjd, "RMSE", "IFN + SJD")
p4 <- plot_knn_metric(knn_ifn_sjd, "MAPE", "IFN + SJD")



#-----------------------------

build_long_df <- function(data, dataset_name) {
  data$folds %>%
    mutate(dataset = dataset_name) %>%
    pivot_longer(
      cols = c(RMSE, MAE),
      names_to = "metric",
      values_to = "value"
    )
}

df_all <- bind_rows(
  build_long_df(knn_ifn, "CTRL + SjD"),
  build_long_df(knn_ifn_sjd, "SjD only")
)

# ordre couleur global
ordre_couleur <- unique(df_all$method)

df_all <- data.table(df_all %>%
  group_by(dataset, metric, method) %>%
  mutate(med = median(value, na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(dataset, metric) %>%
  mutate(
    method = factor(method, levels = unique(method[order(med)]))
  ))


ordre_method_global <- knn_ifn_sjd$folds %>%
  group_by(method) %>%
  summarise(med_rmse = median(RMSE, na.rm = TRUE)) %>%
  arrange(med_rmse) %>%
  pull(method)

df_all <- df_all %>%
  mutate(
    method = factor(method, levels = ordre_method_global)
  )


df_all[,metric:=factor(metric, levels=c("RMSE","MAE"))]

plot=ggplot(df_all, aes(method, value, color = method)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  scale_colour_viridis_d(
    option = "viridis",
    limits = ordre_couleur
  ) +
  facet_grid(dataset ~ metric, scales = "free_y") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "right",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9),
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.6
    ))+
      labs(x=NULL, y=NULL)+
  guides(color = guide_legend(title = "Module Scoring Method"))

  
# ggsave(plot=plot, filename=paste0(folder_for_res, "rmse_mae_ifn.pdf"), width=12, height = 10)


setnames(x = meta_data, old="diag", new="DIAGNOSIS")
hist_ifn=ggplot(meta_data, aes(x = IFN_score, color = DIAGNOSIS, fill = DIAGNOSIS)) +
  geom_density(alpha = 0.3, linewidth = 1) +
  scale_color_manual(
    values = c("Control" = "green", "SjD" = "blue")
  ) +
  scale_fill_manual(
    values = c("Control" = "green", "SjD" = "blue")
  ) +
  theme_minimal()



library(patchwork)

a <- wrap_elements(full = as_grob(plot))
b <- wrap_elements(full = as_grob(hist_ifn))


comb <- a/b  +
  plot_layout(ncol = 2, guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(
    plot.tag          = element_text(size = 20, face = "bold"),
    plot.tag.position = c(0, 1)
  )
ggsave(plot=comb, filename=paste0(folder_for_res, "rmse_mae_ifn.pdf"), width=20, height = 10)

