
#' Prepare folds
#' @param ids patient ids
#' @param k number or folds
#' @param seed seed to set
make_common_folds <- function(ids, k = 10, seed = 123) {
  set.seed(seed)
  n <- length(ids)
  folds <- sample(rep(1:k, length.out = n))
  
  data.table(
    ID = ids,
    fold = folds
  )
}

#' UMAP/PHATE TCGA
#' @param rna_seq_data
#' @param ids patient ids
#' @param pam50 PAM 50 subtypes 
#' @param coldata meta data
#' @param method_name method tested
#' @param n_neighbors number of neighbors
#' @param metric distance metric for umap & phate
make_umap_phate_plots_tcga <- function(
    rna_seq_data,
    ids,
    pam50,
    coldata,
    method_name,
    n_neighbors = 15,
    metric
) {
  
  ## ------------------ UMAP
  umap_dt <- data.table(
    uwot::umap(
      X = rna_seq_data,
      n_neighbors = n_neighbors,
      n_components = 2,
      metric = metric
    )
  )
  setnames(umap_dt, c("UMAP1", "UMAP2"))
  umap_dt[, ID := ids]
  
  umap_dt[pam50, PAM50_subtype := i.BRCA_Subtype_PAM50, on = "ID"]
  umap_dt[coldata, Tissue_status := i.definition, on = "ID"]
  
  
  
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
  
  plot_umap_global <- ggplot(
    umap_dt,
    aes(UMAP1, UMAP2, color = PAM50_subtype)
  ) +
    geom_point(alpha = 0.7) +
    scale_color_brewer(palette = "Set1", na.value = "grey80") +
    theme_minimal() +
    theme(
      axis.text = element_blank(),
      axis.title = element_blank()
    ) +
    ggtitle(paste(method_name, "- UMAP (PAM50)")) + theme_cadre
  
  plot_umap_global_status <- ggplot(
    umap_dt,
    aes(UMAP1, UMAP2, color = Tissue_status)
  ) +
    geom_point(alpha = 0.7) +
    scale_color_brewer(palette = "Set2", na.value = "grey80") +
    theme_minimal() +
    theme(
      axis.text = element_blank(),
      axis.title = element_blank()
    ) +
    ggtitle(paste(method_name, "- UMAP (Tissue)")) + theme_cadre
  
  ## ------------------ PHATE
  phate_dt <- data.table(
    phateR::phate(
      data = rna_seq_data,
      ndim = 2,
      knn = n_neighbors,
      knn.dist.method =metric,
      mds.dist.method = metric,
      decay = 20,
      seed = 123,
      verbose = FALSE
    )$embedding
  )
  setnames(phate_dt, c("PHATE1", "PHATE2"))
  phate_dt[, ID := ids]
  
  phate_dt[pam50, PAM50_subtype := i.BRCA_Subtype_PAM50, on = "ID"]
  phate_dt[coldata, Tissue_status := i.definition, on = "ID"]
  
  plot_phate_global <- ggplot(
    phate_dt,
    aes(PHATE1, PHATE2, color = PAM50_subtype)
  ) +
    geom_point(alpha = 0.7) +
    scale_color_brewer(palette = "Set1", na.value = "grey80") +
    theme_minimal() +
    theme(
      axis.text = element_blank(),
      axis.title = element_blank()
    ) +
    ggtitle(paste(method_name, "- PHATE (PAM50)"))
  
  plot_phate_global_status <- ggplot(
    phate_dt,
    aes(PHATE1, PHATE2, color = Tissue_status)
  ) +
    geom_point(alpha = 0.7) +
    scale_color_brewer(palette = "Set2", na.value = "grey80") +
    theme_minimal() +
    theme(
      axis.text = element_blank(),
      axis.title = element_blank()
    ) +
    ggtitle(paste(method_name, "- PHATE (Tissue)"))
  
  ## Retourner les 4 plots dans l’ordre souhaité
  return(list(
    plot_umap_global,
    plot_umap_global_status,
    plot_phate_global,
    plot_phate_global_status
  ))
}

#' UMAP TCGA
#' @param rna_seq_data
#' @param ids patient ids
#' @param pam50 PAM 50 subtypes 
#' @param n_neighbors number of neighbors
#' @param metric distance metric for umap & phate
make_umap_pam50 <- function(
    module_matrix_list,
    ids,
    pam50,
    n_neighbors = 15,
    metric="correlation"
) {
  
  plot_list = list()
  for(method in names(module_matrix_list)){
    message("Method : ", method)
    
  dt = module_matrix_list[[method]]
  dt$ID = NULL
  
  ## ------------------ UMAP
  umap_dt <- data.table(
    uwot::umap(
      X = dt,
      n_neighbors = n_neighbors,
      n_components = 2,
      metric = metric
    )
  )
  setnames(umap_dt, c("UMAP1", "UMAP2"))
  umap_dt[, ID := ids]
  
  umap_dt[pam50, PAM50_subtype := i.BRCA_Subtype_PAM50, on = "ID"]

  
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
  
  plot_umap <- ggplot(
    umap_dt,
    aes(UMAP1, UMAP2, color = PAM50_subtype)
  ) +
    geom_point(alpha = 0.7) +
    scale_color_brewer(palette = "Set1", na.value = "grey80") +
    theme_minimal() +
    theme(
      axis.text = element_blank(),
      axis.title = element_blank()
    )+
    labs(title=method) + theme_cadre
  
  plot_list = list.append(plot_list, plot_umap)
  }
  names(plot_list) = names(module_matrix_list)
  

  return(plot_list)
}


#' Run KNN and RF classification
#' @param module_matrix_list list of module matrices from all methods
#' @param pam50 PAM 50 subtypes 
#' @param folds_dt folds to run
#' @param train_prop proportion of patients for train (between 0 and 1)
#' @param ntree number of tree for random forrest
#' @param k_knn number of neighbors for KNN
#' @param classifier "rf" or "knn"
#' @param seed seed
pam50_rf_knn_classification <- function(
    module_matrix_list,
    pam50,
    folds_dt,
    train_prop = 0.7,
    ntree = 500,
    k_knn = 20,
    classifier = c("rf", "knn"),
    seed = 123
) {
  
  library(data.table)
  library(randomForest)
  library(pROC)
  library(ggplot2)
  library(caret)
  
  classifier <- match.arg(classifier)
  
  ## ---------------- Sanity checks
  stopifnot("ID" %in% names(pam50))
  stopifnot("BRCA_Subtype_PAM50" %in% names(pam50))
  
  pam50 <- pam50[
    BRCA_Subtype_PAM50 %in% c("Basal", "Her2", "LumA", "LumB")
  ]
  pam50[, BRCA_Subtype_PAM50 := factor(
    BRCA_Subtype_PAM50,
    levels = c("Basal", "Her2", "LumA", "LumB")
  )]
  
  ## ---------------- Containers
  roc_store <- list()
  auc_list  <- list()
  acc_list  <- list()
  
  #--  METHODS LOOP
  for (method_name in names(module_matrix_list)) {
    
    message("Processing method: ", method_name,
            " | classifier: ", classifier)
    
    X <- copy(module_matrix_list[[method_name]])
    
    ## Merge PAM50 + folds
    dt <- merge(X, pam50, by = "ID")
    dt <- merge(dt, folds_dt, by = "ID")
    dt <- dt[!BRCA_Subtype_PAM50 %in% c("Normal", "NA")]
    
    features <- setdiff(
      names(dt),
      c("ID", "BRCA_Subtype_PAM50", "fold")
    )
    
    ## ---------- ACCURACY (multi-class, fold-wise)
    acc_folds <- numeric(length(unique(dt$fold)))
    
    ## ---------- AUC (one-vs-rest)
    auc_by_class <- list()
    
    for (f in sort(unique(dt$fold))) {
      
      set.seed(seed + f)
      
      fold_data <- dt[fold == f]
      
      idx_train <- sample(
        seq_len(nrow(fold_data)),
        size = floor(train_prop * nrow(fold_data))
      )
      
      train <- fold_data[idx_train]
      test  <- fold_data[-idx_train]
      
      ## ----- Factor response
      train[, BRCA_Subtype_PAM50 := factor(
        BRCA_Subtype_PAM50,
        levels = levels(pam50$BRCA_Subtype_PAM50)
      )]
      test[, BRCA_Subtype_PAM50 := factor(
        BRCA_Subtype_PAM50,
        levels = levels(train$BRCA_Subtype_PAM50)
      )]
      
      ## ----- Feature matrices
      x_train <- as.matrix(train[, ..features])
      x_test  <- as.matrix(test[, ..features])
      
      ## --MULTI-CLASS ACCURACY--
      if (classifier == "rf") {
        
        model_mc <- randomForest(
          x = x_train,
          y = train$BRCA_Subtype_PAM50,
          ntree = ntree
        )
        
        pred_mc <- predict(model_mc, x_test, type = "class")
        
      } else {
        
        model_mc <- knn3(
          x = x_train,
          y = train$BRCA_Subtype_PAM50,
          k = k_knn
        )
        
        pred_mc <- predict(model_mc, x_test, type = "class")
      }
      
      acc_folds[f] <- mean(pred_mc == test$BRCA_Subtype_PAM50)
      
      #--- ONE-vs-REST AUC--
      for (cls in levels(train$BRCA_Subtype_PAM50)) {
        
        y_train_bin <- factor(train$BRCA_Subtype_PAM50 == cls)
        y_test_bin  <- factor(test$BRCA_Subtype_PAM50 == cls)
        
        if (classifier == "rf") {
          
          model_bin <- randomForest(
            x = x_train,
            y = y_train_bin,
            ntree = ntree
          )
          
          prob <- predict(
            model_bin,
            x_test,
            type = "prob"
          )[,"TRUE"]
          
        } else {
          
          model_bin <- knn3(
            x = x_train,
            y = y_train_bin,
            k = k_knn
          )
          
          prob <- predict(
            model_bin,
            x_test,
            type = "prob"
          )[,"TRUE"]
        }
        
        if (is.null(auc_by_class[[cls]])) {
          auc_by_class[[cls]] <- list(y = c(), prob = c())
        }
        
        auc_by_class[[cls]]$y    <- c(
          auc_by_class[[cls]]$y,
          as.numeric(y_test_bin) - 1
        )
        auc_by_class[[cls]]$prob <- c(
          auc_by_class[[cls]]$prob,
          prob
        )
      }
    }
    
    ## ---------- Final ACCURACY + CI
    acc_mean <- mean(acc_folds)
    acc_sd   <- sd(acc_folds)
    n_folds  <- length(acc_folds)
    
    acc_ci_low  <- acc_mean - 1.96 * acc_sd / sqrt(n_folds)
    acc_ci_high <- acc_mean + 1.96 * acc_sd / sqrt(n_folds)
    
    acc_list[[method_name]] <- data.table(
      Method = method_name,
      Accuracy = acc_mean,
      Accuracy_CI_low = acc_ci_low,
      Accuracy_CI_high = acc_ci_high
    )
    
    ## ---------- Final AUC + ROC storage
    for (cls in names(auc_by_class)) {
      
      roc_obj <- roc(
        auc_by_class[[cls]]$y,
        auc_by_class[[cls]]$prob,
        quiet = TRUE
      )
      
      auc_ci <- ci.auc(roc_obj)
      
      auc_list[[paste(method_name, cls, classifier, sep = "_")]] <-
        data.table(
          Method = method_name,
          Class  = cls,
          Classifier = classifier,
          AUC    = as.numeric(auc(roc_obj)),
          CI_low = auc_ci[1],
          CI_high= auc_ci[3]
        )
      
      roc_store[[cls]][[paste(method_name, classifier, sep = "_")]] <- roc_obj
    }
  }
  
  ##--- ROC PLOTS (1 per class) ---
  roc_plots_by_class <- list()
  
  for (cls in names(roc_store)) {
    
    roc_dt <- rbindlist(
      lapply(names(roc_store[[cls]]), function(key) {
        
        roc_obj <- roc_store[[cls]][[key]]
        
        data.table(
          Class = cls,
          Method = key,
          FPR = 1 - roc_obj$specificities,
          TPR = roc_obj$sensitivities
        )
      })
    )
    
    method_order <- c(
      "GEDO_knn",
      "GEDOcorr_knn",
      "UMAP_GEDO_knn",
      "PCA1_knn",
      "MZS_knn",
      "ssGSEA_knn",
      "GSVA_knn"
    )
    roc_dt[, Method := factor(Method, levels = method_order)]
    roc_dt[,Method:=gsub(x = Method, pattern = "_knn", replacement = "")]
    roc_dt[,Method:=factor(Method, levels=c("GEDO","GEDOcorr","UMAP_GEDO","PCA1","MZS","ssGSEA","GSVA"))]
    
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
    
    roc_plots_by_class[[cls]] <-
      ggplot(roc_dt, aes(FPR, TPR, color = Method)) +
      geom_line(linewidth = 1) +
      geom_abline(linetype = "dashed", color = "grey50") +
      scale_colour_viridis_d(
        option = "viridis",
        name = "Module Scoring Method"
      )+
      theme_minimal() +
      labs(
        title = cls,
        x = "False Positive Rate",
        y = "True Positive Rate"
      )+
      coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, clip = "on") + theme_cadre
    
  }
  
  
  # FORMAT TABLE (AUC + ACCURACY)
  
  auc_dt <- rbindlist(auc_list)
  
  auc_dt[, Class := factor(Class,
                           levels = c("LumA", "LumB", "Her2", "Basal"))]
  auc_dt[, Method := factor(Method,
                            levels = c("GEDO","GEDOcorr","UMAP_GEDO",
                                       "PCA1","MZS","ssGSEA","GSVA"))]
  
  ## AUC formatting
  auc_dt[, AUC_fmt := formatC(AUC, digits = 3, format = "fg")]
  auc_dt[, CI_low_fmt := formatC(CI_low, digits = 3, format = "fg")]
  auc_dt[, CI_high_fmt := formatC(CI_high, digits = 3, format = "fg")]
  
  auc_dt[, AUC_CI := paste0(
    AUC_fmt, " [", CI_low_fmt, " ; ", CI_high_fmt, "]"
  )]
  
  auc_wide <- dcast(
    auc_dt,
    Method ~ Class,
    value.var = "AUC_CI"
  )
  
  ## Accuracy formatting
  acc_dt <- rbindlist(acc_list)
  
  acc_dt[, Method := factor(Method,
                            levels = c("GEDO","GEDOcorr","UMAP_GEDO",
                                       "PCA1","MZS","ssGSEA","GSVA"))]
  
  acc_dt[, Acc_fmt := formatC(Accuracy, digits = 3, format = "fg")]
  acc_dt[, CI_low_fmt := formatC(Accuracy_CI_low, digits = 3, format = "fg")]
  acc_dt[, CI_high_fmt := formatC(Accuracy_CI_high, digits = 3, format = "fg")]
  
  acc_dt[, Accuracy_CI := paste0(
    Acc_fmt, " [", CI_low_fmt, " ; ", CI_high_fmt, "]"
  )]
  
  acc_final <- acc_dt[, .(Method, Accuracy = Accuracy_CI)]
  
  ## Merge AUC + Accuracy
  table_article <- merge(
    auc_wide,
    acc_final,
    by = "Method"
  )
  
  ## ---- RETURN ---
  list(
    roc_plots_by_class = roc_plots_by_class,
    auc_table          = auc_dt,
    accuracy_table     = acc_dt,
    formatted_table    = table_article,
    roc_store          = roc_store
  )
}


#' Run ordinarity metrics across PAM50 subtypes 
#' @param module_matrix_list list of module matrices from all methods
#' @param pam50 PAM 50 subtypes 
#' @param knn number of neighbors 
#' @param order_levels vector with the levels of categories to test
#' @param seed seed
evaluate_pam50_ordinality_phate <- function(
    module_matrix_list,
    pam50,
    order_levels = c("LumA", "LumB", "Her2", "Basal"),
    knn = 15,
    seed = 123
) {
  

  
  ## ---------- Ordre biologique
  order_code <- setNames(seq_along(order_levels), order_levels)
  
  ## ---------- Containers
  results_list <- list()
  plot_list <- list()
  
  #METHODS LOOP
  for (method_name in names(module_matrix_list)) {
    
    message("Processing ordinality analysis for method: ", method_name)
    
    X <- copy(module_matrix_list[[method_name]])
    
    ## ---------- Merge PAM50
    dt <- merge(X, pam50, by = "ID")
    dt <- dt[BRCA_Subtype_PAM50 %in% order_levels]
    
    dt[, subtype_ord := order_code[BRCA_Subtype_PAM50]]
    
    ## ---------- PHATE
    # message("Processing ordinality analysis for method: ", method_name, " : PHATE...")
    
    set.seed(seed)
    ph <- phate(
      data = as.matrix(dt[, !c("ID", "BRCA_Subtype_PAM50", "subtype_ord"), with = FALSE]),
      ndim = 1,
      knn = knn,
      verbose = FALSE
    )
    
    dt[, PHATE1 := ph$embedding[, 1]]
    med_luma  <- median(dt[BRCA_Subtype_PAM50 == "LumA", PHATE1])
    med_basal <- median(dt[BRCA_Subtype_PAM50 == "Basal", PHATE1])
    
    if (med_luma > med_basal) {
      dt[, PHATE1 := -PHATE1]
    }
    
    
    ## ---------- Kendall tau
    tau <- cor(
      dt$PHATE1,
      dt$subtype_ord,
      method = "kendall"
    )
    
    ## ---------- Jonckheere–Terpstra
    dt[, BRCA_Subtype_PAM50:=factor(BRCA_Subtype_PAM50, levels = order_levels)]
    dt[,BRCA_Subtype_PAM50_num := as.numeric(BRCA_Subtype_PAM50)]
    jt <- jonckheere.test(x = dt$PHATE1, g=dt$BRCA_Subtype_PAM50_num, alternative = "increasing")
    
   
    
    ## ---------- Store results
    results_list[[method_name]] <- data.table(
      Method = method_name,
      Kendall_tau = tau,
      JT_statistic = unname(jt$statistic),
      JT_pvalue = jt$p.value
    )
    
    ## ---------- Plot
    pval_txt <- ifelse(
      jt$p.value < 1e-16,
      "< 1e-16",
      formatC(jt$p.value, format = "e", digits = 2)
    )
    
    subtitle_txt <- sprintf(
      "Kendall tau = %.2f | JT p = %s",
      tau,
      pval_txt
    )
    
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
    
    p <- ggplot(
      dt,
      aes(
        x = factor(BRCA_Subtype_PAM50, levels = order_levels),
        y = PHATE1, 
        color=BRCA_Subtype_PAM50
    )) +
      geom_boxplot(outlier.shape = NA, alpha = 0.7) +
      geom_jitter(width = 0.2, size = 0.6, alpha = 0.5) +
      scale_color_brewer(palette = "Set1", na.value = "grey80") +
      theme_minimal() +
      theme(
        legend.position = "none"
      ) +
      labs(
        title = method_name,
        subtitle = subtitle_txt,
        x = "PAM50 subtype (ordered)",
        y = "PHATE1"
      ) + theme_cadre
    
    plot_list[[method_name]] <- p
  }
  
  ## ---------- Final table
  results_table <- rbindlist(results_list)
  
  ## ---------- Combined figure
  combined_plot <- wrap_plots(plot_list, ncol = 2)
  
  ## ---------- Return
  return(list(
    results_table = results_table,
    plots_by_method = plot_list,
    combined_plot = combined_plot
  ))
}


