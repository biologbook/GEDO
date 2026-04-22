
#' load packages and install if needed
#' @param packages_list character vector of packages names to load
install_and_load <- function(packages_list) {
  not_installed <- packages_list[!packages_list %in% installed.packages()[, "Package"]]
  not_installed_bioconductor=not_installed[not_installed %in% c("msigdbr")]
  not_installed_other = not_installed[!not_installed %in% not_installed_bioconductor] 
  
  if (length(not_installed_other) > 0) {
    install.packages(not_installed)
  }
  
  if (length(not_installed_bioconductor) > 0) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install(not_installed_bioconductor)
  }
  if("msigdbdf" %in% not_installed){
    install.packages('msigdbdf', repos = c('https://igordot.r-universe.dev', 'https://cloud.r-project.org'))
  }
  print(sapply(packages_list, require, character.only = TRUE))
}


detach_all_packages <- function() {
  basic_packages <- c("package:stats", "package:graphics", "package:grDevices",
                      "package:utils", "package:datasets", "package:methods", "package:base")
  
  loaded_packages <- search()[grepl("^package:", search())]
  to_detach <- setdiff(loaded_packages, basic_packages)
  
  for (pkg in to_detach) {
    detach(pkg, character.only = TRUE, unload = TRUE)
  }
}



heatmap.ifn=function(gedo_obj, title="", IFN_score, legend=T, annotation_legend=T){
  #1. preparation of data
  diag = gedo_obj$diag
  module_matrix = as.matrix(gedo_obj$module_matrix)
  id=1:nrow(module_matrix)
  module_matrix=t(module_matrix)
  colnames(module_matrix)=id

  #2. Legends


  #
  # IFN_score[!is.finite(IFN_score)] = NA
  # annotation_col=data.frame(diag=diag, IFN_score=IFN_score)
  # rownames(annotation_col)=id
  #
  reference_group=as.character(gedo_obj$config$reference_group)
  another_group = as.character(unique(diag[!diag %in% gedo_obj$config$reference_group]))
  #
  # overlap_colors <- list(
  #   diag = setNames(c("green", "blue"), c(reference_group, another_group)),
  #   IFN_score=c("blue", "white", "red")
  # )


  IFN_score[!is.finite(IFN_score)] <- NA

  annotation_col <- data.frame(diag = diag, IFN_score = IFN_score)
  rownames(annotation_col) <- colnames(module_matrix)

  # 2) Ne PAS donner de couleurs pour IFN_score (numérique)
  overlap_colors <- list(
    diag = setNames(c("green", "blue"),
                    c(reference_group, another_group))
  )

  # IFN_score : ne pas donner un simple vecteur non nommé,
  # mais un dégradé généré par colorRampPalette
  # overlap_colors$IFN_score <- colorRampPalette(c("blue", "white", "red"))(100)



  #3. Heatmap
  heatmap=pheatmap(
    module_matrix,
    annotation_col = annotation_col,
    annotation_colors = overlap_colors,
    cluster_rows = T,
    cluster_cols = T,
    show_rownames = FALSE,
    show_colnames = FALSE,
    scale = "none",
    annotation_names_row=F,
    annotation_legend=annotation_legend,
    legend=legend
  )
  return(heatmap)
}
# 
# 
# heatmap.ifn_XY=function(gedo_obj, title="", IFN_score, legend=T, annotation_legend=T){
#   #1. preparation of data
#   diag = gedo_obj$diag
#   module_matrix = as.matrix(gedo_obj$module_matrix)
#   cols= colnames(module_matrix)
#   
#   id=1:nrow(module_matrix)
#   module_matrix=t(module_matrix)
#   colnames(module_matrix)=id
#   
#   
# omic = ifelse(test = grepl("^X", cols), yes = "T","M")
#   
#   #2. Legends
#   
#   
#   # 
#   # IFN_score[!is.finite(IFN_score)] = NA
#   # annotation_col=data.frame(diag=diag, IFN_score=IFN_score)
#   # rownames(annotation_col)=id
#   # 
#   reference_group=as.character(gedo_obj$config$reference_group)
#   another_group = as.character(unique(diag[!diag %in% gedo_obj$config$reference_group]))
#   # 
#   # overlap_colors <- list(
#   #   diag = setNames(c("green", "blue"), c(reference_group, another_group)),
#   #   IFN_score=c("blue", "white", "red")
#   # )
#   
#   
#   IFN_score[!is.finite(IFN_score)] <- NA
#   
#   annotation_col <- data.frame(diag = diag, IFN_score = IFN_score)
#   rownames(annotation_col) <- colnames(module_matrix)
#   
#   annotation_row = data.frame(omic = omic)
#   rownames(annotation_row)=cols
#   
#   # 2) Ne PAS donner de couleurs pour IFN_score (numérique)
#   overlap_colors <- list(
#     diag = setNames(c("green", "blue"),
#                     c(reference_group, another_group))
#   )
#   
#   # IFN_score : ne pas donner un simple vecteur non nommé,
#   # mais un dégradé généré par colorRampPalette
#   # overlap_colors$IFN_score <- colorRampPalette(c("blue", "white", "red"))(100)
#   
#   
#   
#   #3. Heatmap
#   heatmap=pheatmap(
#     module_matrix,
#     annotation_col = annotation_col, 
#     annotation_row=annotation_row,
#     annotation_colors = overlap_colors,
#     cluster_rows = T,               
#     cluster_cols = T,               
#     show_rownames = FALSE,              
#     show_colnames = FALSE,             
#     scale = "none",                    
#     annotation_names_row=F,
#     annotation_legend=annotation_legend,
#     legend=legend
#   )
#   return(heatmap)
# }
# 



#' compute PCA1 on one gene module
#' @param data data.table of gene matrix for one module (M). Genes in columns, individuals in line. Only Gene columns, no id. 
#' @param module_name gene module in GSEA database (extraction with msigdbr, column gs_name).
#' @param diag character vector with two groups : Controls and Diseased. diag must correspond to lines in data. 
#' @param reference_group character with the name of control group in diag.

compute_pca1=function(data, module_name, diag, reference_group, category, subcategory, charge_reactome_modules=F){ 
  
  
  #0. if necessary, chargind reactome modules
  if(charge_reactome_modules==T){
    reactome_modules <- data.table(msigdbr(species = "Homo sapiens", collection = category, subcollection = subcategory))
    reactome_modules <- reactome_modules[, .SD, .SDcols = c("gs_name", "ensembl_gene")]
    module_number = length(unique(reactome_modules$gs_name))
    modules = unique(reactome_modules$gs_name)
  }               
  
  #1. Selecting data for the module
  genes_in_module =  unique(reactome_modules[gs_name == module_name]$ensembl_gene)
  genes_in_data = colnames(data)
  genes_to_keep = genes_in_data[genes_in_data %in% genes_in_module]
  if(length(genes_to_keep)<2){stop("no genes in data for this module")}
  data_module = data[,..genes_to_keep, with=F]
  
  #Delete columns with only one value
  cols_to_keep <- names(data_module)[sapply(data_module, function(x) uniqueN(x, na.rm = TRUE) > 1)]
  data_module <- data_module[, ..cols_to_keep]
  
  
  
  #2. PCA
  pca_result <- prcomp(data_module, center = TRUE, scale. = TRUE)
  pc1_values <- pca_result$x[, 1]
  
  table = data.table(pca1_values=pc1_values, diag=diag)
  mean_pca1_groupe1=mean(table[diag==reference_group]$pca1_values)
  mean_pca1_groupe2=mean(table[diag!=reference_group]$pca1_values)
  
  
  #3. If PCA1 for Control is superior of PCA1 for Diseased the sign is inverted.
  if(mean_pca1_groupe1 > mean_pca1_groupe2){
    pc1_values <- -pc1_values
  }
  
  pc1_values = (pc1_values - min(pc1_values)) / (max(pc1_values) - min(pc1_values))
  
  return(pc1_values)
}


#' compute mean of Z-scores on one gene module
#' @param data data.table of gene matrix for one module (M). Genes in columns, individuals in line. Only Gene columns, no id. 
#' @param module_name gene module in GSEA database (extraction with msigdbr, column gs_name).
#' @param diag character vector with two groups : Controls and Diseased. diag must correspond to lines in data. 
#' @param reference_group character with the name of control group in diag.
compute_mean_z_score=function(data, module_name, diag, reference_group){ 
  
  
  #1. Selecting data for the module
  genes_in_module =  unique(reactome_modules[gs_name == module_name]$ensembl_gene)
  genes_in_data = colnames(data)
  genes_to_keep = genes_in_data[genes_in_data %in% genes_in_module]
  if(length(genes_to_keep)<2){stop()}
  data_module = data[,..genes_to_keep, with=F]
  
  #2. Mean and SD of controls
  data_module$diag = diag
  
  gr1 = data_module[diag==reference_group]
  gr1[,diag:=NULL]
  data_module[,diag:=NULL]
  
  gr1_means = colMeans(gr1)
  gr1_sds = apply(gr1, 2, sd)
  
  #3. Z scores
  # calculate_z_scores <- function(sample, means, sds) {
  #   (sample - means) / sds
  # }
  calculate_z_scores <- function(sample, means, sds) {
    z <- numeric(length(sample))
    
    zero_sd <- sds == 0 | is.na(sds)
    
    # cas normal
    z[!zero_sd] <- (sample[!zero_sd] - means[!zero_sd]) / sds[!zero_sd]
    
    # cas sd == 0
    z[zero_sd] <- means[zero_sd]
    
    z
  }
  
  
  
  z_scores <- apply(data_module, 1, calculate_z_scores, means = gr1_means, sds = gr1_sds)
  
  #Force result in matrix
  if (is.null(dim(z_scores))) {
    z_scores <- matrix(z_scores, ncol = 1)
  } else {
    z_scores <- t(z_scores)
  }
  
  
  mean_z_scores <- rowMeans(z_scores, na.rm = TRUE)
  
  mean_z_scores = (mean_z_scores - min(mean_z_scores)) / (max(mean_z_scores) - min(mean_z_scores))
  
  return(mean_z_scores)
}






#' Plot AUC of each module from GEDO object
#' @param gedo_obj GEDO object computed with gedo()
gedo.module.auc=function(gedo_obj){
  module_matrix=gedo_obj$module_matrix
  diag = gedo_obj$diag
  
  #AUC table
  auc_results <- rbindlist(lapply(colnames(module_matrix), function(module) {
    auc_value <- pROC::auc(diag, module_matrix[[module]])
    
    data.table(
      gene_module = module,
      auc = as.numeric(auc_value)
    )
  }))
  auc_results=auc_results[order(-auc)]
  
  #AUc plot
  
  min_auc_top_50 = min(auc_results[1:50]$auc)
  plot_50=ggplot(auc_results[1:50], aes(x = reorder(gene_module, auc))) +
    geom_segment(aes(xend = reorder(gene_module, auc),
                     y = min_auc_top_50, yend = auc),
                 color = "steelblue", linewidth = 3) +
    coord_flip() +
    theme_minimal() +
    scale_y_continuous(limits = c(min_auc_top_50, 1)) +
    labs(title="50 modules with the best AUC",
         x = "Gene modules",
         y = "AUC"
    ) +
    theme(
      axis.text.y = element_text(size = 8),
      plot.title = element_text(hjust = 0.5)
    )
  
  
  
  min_auc = min(auc_results$auc)
  
  plot_all=ggplot(auc_results, aes(x = reorder(gene_module, auc))) +
    geom_segment(aes(
      xend = reorder(gene_module, auc),
      y = min_auc,
      yend = auc
    ),
    color = "steelblue", linewidth = 3) +
    coord_flip() +
    scale_y_continuous(limits = c(min_auc, 1)) +
    labs(title = "All Modules",
         x = "Gene modules",
         y = "AUC"
    ) +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.major.y = element_line(color = "grey80"),
      plot.title = element_text(hjust = 0.5)
    )
  
  
  combined_plot= plot_all + plot_50 + plot_layout(widths = c(1, 2))
  print(combined_plot)
  
  list_to_return=list(auc_results=auc_results, plot_auc=combined_plot)
  
  return(list_to_return)
}








#' compute Module matrix with PCA1 or mean of z-scores from RNAseq data
#' @param method method for scoring gene modules : 'mean_z_score' or 'pca1'
#' @param data data.table of gene matrix for one module (M). Genes in columns, individuals in line. Only Gene columns, no id. 
#' @param diag character vector with two groups : Controls and Diseased. diag must correspond to lines in data. 
#' @param reference_group character with the name of control group in diag.
#' @param category if  reactome_modules not provided, the collection in GSEA for extraction with with msigdbr
#' @param subcategory if reactome_modules not provided, the subcollection in GSEA for extraction with with msigdbr
#' @param num_cores Number of cores for parallelisation. By default, num_cores <- detectCores() - 1. 

compute_module_matrix = function(method, data, diag, reference_group, category, subcategory, num_cores=NULL){
  if(is.null(num_cores)){num_cores <- detectCores() - 1}
  
  packages_list = c("data.table", "FNN", "magrittr", "tryCatchLog","rlist",
                    "igraph", "pbapply","rgl","rdist","msigdbr","msigdbdf",
                    "future.apply","dbscan", "progressr",
                    "gridExtra","MASS","ggplot2", "patchwork", "dplyr", "tidyr","parallel","reticulate", "tryCatchLog")
  #1. Config file
  config=list(
    groups= table(diag),
    reference_group=reference_group,
    category=category, 
    subcategory=subcategory, 
    num_cores=num_cores,
    method=method
  )
  
  #2. Loading gene modules
  reactome_modules <- data.table(msigdbr(species = "Homo sapiens", collection = category, subcollection  = subcategory))
  reactome_modules <- reactome_modules[, .SD, .SDcols = c("gs_name", "ensembl_gene")]
  module_number = length(unique(reactome_modules$gs_name))
  modules = unique(reactome_modules$gs_name)
  
  
  #3. Verifications on diag vector
  if(length(unique(diag))>2){
    print("Error : more than two classes in diag") 
    return(NULL)}
  if(length(unique(diag))<2){
    print("Error : less than two classes in diag") 
    return(NULL)}
  if(any(is.na(diag))){
    print("Error : NA in diag") 
    return(NULL)}
  
  
  #4. Computing PCA1 or mean of z-scores for each modules
  
  plan(multicore, workers = num_cores)
  handlers("txtprogressbar") 
  cat("Parallelization on modules","\n")
  
  export=list(modules=modules,reactome_modules=reactome_modules,method=method, data=data, diag=diag, reference_group=reference_group, category=category, subcategory=subcategory, compute_pca1=compute_pca1, compute_mean_z_score=compute_mean_z_score)
  
  #4.1. PCA1
  if(method=="pca1"){
    module_matrix_list = with_progress({
      p <- progressor(along = modules)
      
      future_lapply(modules, function(module_name) {
        p(message = sprintf("Module: %s", module_name))
        library(tryCatchLog)
        
        
        tryCatch({
          res = as.numeric(compute_pca1(
            data = data,
            module_name = module_name,
            diag = diag,
            reference_group = reference_group))
          return(res)
        }, error = function(e) {
          cat("Error in module:", module_name, "\nMessage:", e$message, "\n")
          return(NULL)
        })
      }, future.seed=T, future.globals = c(export, list(p=p)), future.packages = packages_list)
    })
  }
  
  
  #4.2. mean of z-scores : 
  if(method=="mean_z_score"){
    module_matrix_list = with_progress({
      p <- progressor(along = modules)
      
      future_lapply(modules, function(module_name) {
        p(message = sprintf("Module: %s", module_name))
        library(tryCatchLog)
        
        
        tryCatch({
          res = as.numeric(compute_mean_z_score(
            data = data,
            module_name = module_name,
            diag = diag,
            reference_group = reference_group))
          
          return(res)
        }, error = function(e) {
          cat("Error in module:", module_name, "\nMessage:", e$message, "\n")
          return(NULL)
        })
      }, future.seed=T, future.globals = c(export, list(p=p)), future.packages = packages_list)
    })
  }
  
  names(module_matrix_list)=modules
  module_matrix_list <- Filter(Negate(is.null), module_matrix_list)
  module_matrix <- data.table(do.call(cbind, module_matrix_list))
  # cols_with_na <- sapply(module_matrix, function(col) any(is.na(unlist(col))))
  # module_matrix[, names(module_matrix)[cols_with_na] := NULL]
  
  
  
  
  #5. Return gedo-like object
  module_matrix_obj=list(
    module_matrix = module_matrix,
    config=config,
    diag=diag
  )
  return(module_matrix_obj)
}




#' compute comparison of module AUC to predict Control/Diseased, according to each method (GEDO, PCA1, MEAN OF Z-SCORES)
#' @param matrix_list list with objects from gedo() or compute_module_matrix(). the list must be named.
compute_auc_modules = function(matrix_list){
  
  #1. Computing AUC for each module and each element of matrix_list
  auc_results <- list()
  
  for(metric_name in names(matrix_list)){
    print(metric_name)
    obj = matrix_list[[metric_name]]
    dt = obj$module_matrix
    modules = colnames(dt)
    dt$diag = obj$diag
    
    for (pathway in setdiff(x = colnames(dt), y = c("diag"))) {
      print(paste("method=", metric_name, ". module=", pathway))
      auc_value <- auc(dt$diag, dt[[pathway]])
      
      auc_results <- append(auc_results, list(data.table(
        metric = metric_name,
        pathway = pathway,
        auc = as.numeric(auc_value)
      )))
    }
    
  }
  dt_auc_results <- rbindlist(auc_results, use.names = TRUE, fill = TRUE)
  
  dt_auc_results[, median_auc := median(auc), by=metric]
  dt_auc_results <- dt_auc_results[order(-median_auc)]
  best_metric=unique(dt_auc_results[median_auc==max(median_auc)]$metric)
  metrics_other_the_best = unique(dt_auc_results[metric!=best_metric]$metric)
  
  levels=rev(unique(dt_auc_results$metric))
  dt_auc_results[, metric := factor(metric, levels = levels)]
  
  dt_auc_results[, group := ifelse(metric == best_metric, best_metric, "Others")]
  dt_auc_results[, is_gedo := grepl("gedo", metric, ignore.case = TRUE)]
  
  comparisons_list <- lapply(metrics_other_the_best, function(m) c(best_metric, m))
  
  
  #2. Boxplot
  boxplot = ggplot(dt_auc_results, aes(x = metric, y = auc, color = is_gedo)) +
    geom_boxplot() +
    geom_jitter(alpha = 0.05) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_color_manual(values = c("TRUE" = "firebrick", "FALSE" = "black")) +
    guides(color = "none") +  
    stat_compare_means(
      comparisons = comparisons_list,
      method = "wilcox.test",
      paired = TRUE,
      label = "p.signif",
      hide.ns = FALSE,
      p.adjust.method = "BH"   # ou "bonferroni", "holm", etc.
    ) +
    labs(x="Module Scoring Method", y="AUC", title="a")+
    coord_flip()
  
  
  #3. AUC face to face
  auc_wide <- dcast(dt_auc_results, pathway ~ metric, value.var = "auc")
  
  plot_list=list()
  
  for(m in metrics_other_the_best){
    
    plot=ggplot(auc_wide, aes(x = .data[[m]], y = .data[[best_metric]], label = pathway)) +
      geom_point(size = 1, color = "black", alpha=0.3) +  # Points pour chaque pathway
      theme_minimal() +
      geom_abline(slope = 1, intercept = 0, color = "red", size = 1)+
      geom_quantile(quantiles = 0.5, method = "rqss", lambda = 0.1, color = "blue", size = 1)+  
      xlim(c(0.5,1))+
      ylim(c(0.5,1))+
      coord_fixed() 
    if(m ==metrics_other_the_best[[1]]){
      plot=plot+labs(
        title = "b",
        x = paste0("AUC (",m,")"),
        y = paste0("AUC (", best_metric, ")")
      )}else{
        plot=plot+labs(
          x = paste0("AUC (",m,")"),
          y = paste0("AUC (", best_metric, ")")
        ) 
      }
    
    
    plot_list=list.append(plot_list, plot) 
  }
  
  #4. return combined plot
  combined_subplots <- Reduce(`+`, plot_list) + plot_layout(ncol = length(metrics_other_the_best)/2)
  combined_plot <- boxplot / (Reduce(`+`, plot_list) + plot_layout(ncol = length(matrix_list)))
  list_to_return=list(dt_auc_results=dt_auc_results, auc_wide=auc_wide, plot=combined_plot, boxplot=boxplot, combined_subplots=combined_subplots)
  
  return(list_to_return)
}



#' compute prediction performance of each Module Matrices with Random Forest or KNN
#' @param model The model for prediction : "knn" or "rf"
#' @param k the number of neighbors for knn prediction. 
#' @param k_folds the number of folds for cross validation 
#' @param dt_list list with objects from gedo() or compute_module_matrix(). the list must be named.
#' @param num_cores Number of cores for parallelisation of folds. By default, num_cores <- detectCores() - 1. 

# compute_prediction <- function(model = "knn", k, k_folds = 10, dt_list, num_cores=NULL) {
#   if(is.null(num_cores)){num_cores <- detectCores() - 1}
#   
#   
#   #1. computing k-folds
#   print("computing k-folds")
#   roc_curves <- lapply(names(dt_list), function(dt_name) {
#     cat(paste0(dt_name, "\n"))
#     dt <- dt_list[[dt_name]]
#     diag = dt$diag
#     dt=dt$module_matrix
#     # modules=colnames(dt)
#     dt$diag=diag
#     
#     set.seed(423)
#     trainIndex <- createDataPartition(dt$diag, p = 0.7, list = FALSE)
#     trainData <- dt[trainIndex]
#     testData <- dt[-trainIndex]
#     
#     set.seed(423)
#     folds <- createFolds(trainData$diag, k = k_folds, list = TRUE)
#     
#     all_probs <- c()
#     all_labels <- c()
#     
#     # Parallélisation des k-folds
#     
#     plan(multisession, workers = num_cores)
#     
#     # Progression bar
#     handlers("txtprogressbar") 
#     
#     cat(paste0("Metric : ",dt_name) ,"\n")
#     cat("Parallelization on folds","\n")
#     
#     export = list(model = model, k=k, k_folds =k_folds, folds=folds, trainData=trainData)
# 
#     fold_results = with_progress({
#       p <- progressor(along = folds)
#       
#       future_lapply(1:length(folds), function(i) {
#         cat(paste0("Fold : ",i) ,"\n")
#         
#         testIndexes <- folds[[i]]
#         testFold <- trainData[testIndexes]
#         trainFold <- trainData[-testIndexes]
#         
#         if (model == "knn") {
#           # print("model=knn")
#           classifier <- knn3(diag ~ ., data = trainFold, k = k)
#         } else if (model == "rf") {
#           classifier <- randomForest(as.formula(diag ~ .), data = trainFold, ntree = 400, probability = TRUE)
#         } else {
#           stop("Unsupported model type. Choose 'knn' or 'rf'.")
#         }
#         
#         probs <- predict(classifier, newdata = testFold, type = "prob")[,2]
#         
#         return(list(probs = probs, labels = testFold$diag))
#         
#       }, future.seed = T, future.globals = export, future.packages = c("randomForest","caret"))
#     })
#     
#     
#     
#     all_probs <- unlist(lapply(fold_results, `[[`, "probs"))
#     all_labels <- unlist(lapply(fold_results, `[[`, "labels"))
#     
#     roc_curve <- roc(all_labels, all_probs)
#     cat(paste0("AUC pour ", dt_name, " : ", auc(roc_curve), "\n"))
#     return(roc_curve)
#   })
#   names(roc_curves) <- names(dt_list)
#   
#   print("roc_curves_ok")
# 
#   #2. results in table
#   res_table <- data.table(methods = names(dt_list), AUC = "")
#   for (method in names(roc_curves)) {
#     auc <- roc_curves[[method]]$auc
#     res_table[methods == method, AUC := auc]
#   }
#   setorder(res_table, -AUC)
#   
#   
#   #3. ploting roc curves 
#   print("printing roc curves")
#   roc_data <- do.call(rbind, mclapply(names(roc_curves), function(dt_name) {
#     data.frame(
#       inv_specificity = 1 - roc_curves[[dt_name]]$specificities,
#       sensibility = roc_curves[[dt_name]]$sensitivities,
#       Dataset = dt_name
#     )
#   }, mc.cores = num_cores))
#   
#   roc_data <- data.table(roc_data)
#   desired_order <- names(dt_list)
#   roc_data[, Dataset := factor(Dataset, levels = desired_order)]
#   
#   if(model=="rf"){
#     plot <- ggplot(roc_data, aes(x = inv_specificity, y = sensibility, color = Dataset)) +
#     geom_line(size = 0.5) +
#     labs(title = "(A)", x = "1-specificity", y = "sensibility", color = "Module Scoring Method") +
#     theme_minimal()+
#     guides(color = "none")
#   }
#   if(model=="knn"){
#     plot <- ggplot(roc_data, aes(x = inv_specificity, y = sensibility, color = Dataset)) +
#       geom_line(size = 0.5) +
#       labs(title = "(B)", x = "1-specificity", y = "sensibility", color = "Module Scoring Method") +
#       theme_minimal()
#   }
#   
#   
#   list_to_return = list(auc_results=res_table, roc_curves = plot)
#   
#   return(list_to_return)
# 
# }





compute_prediction_with_ci <- function(model = "knn", k, k_folds = 10, dt_list, num_cores = NULL) {
  # Set default number of cores if not specified
  if (is.null(num_cores)) {
    num_cores <- detectCores() - 1
  }
  
  print("computing k-folds")
  
  roc_curves <- lapply(names(dt_list), function(dt_name) {
    cat(paste0(dt_name, "\n"))
    
    # Load the dataset and add target to the feature matrix
    dt <- dt_list[[dt_name]]
    diag <- dt$diag
    dt <- dt$module_matrix
    dt$diag <- diag
    
    # Split 70% train / 30% test
    set.seed(423)
    trainIndex <- createDataPartition(dt$diag, p = 0.7, list = FALSE)
    trainData <- dt[trainIndex, ]
    testData <- dt[-trainIndex, ]
    
    # Create k folds from the training data
    set.seed(423)
    folds <- createFolds(trainData$diag, k = k_folds, list = TRUE)
    
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
          classifier <- knn3(diag ~ ., data = trainFold, k = k)
        } else if (model == "rf") {
          classifier <- randomForest(as.formula(diag ~ .), data = trainFold, ntree = 400, probability = TRUE)
        } else {
          stop("Unsupported model type. Choose 'knn' or 'rf'.")
        }
        
        # Predict probabilities on test fold
        probs <- predict(classifier, newdata = testFold, type = "prob")[, 2]
        return(list(probs = probs, labels = testFold$diag))
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
      labs(title = "c", x = "1 - Specificity", y = "Sensitivity", color = "Module Scoring Method") +
      theme_minimal() +
      guides(color = "none")
  }
  if (model == "knn") {
    plot <- ggplot(roc_data, aes(x = inv_specificity, y = sensibility, color = Dataset)) +
      geom_line(size = 0.5) +
      scale_colour_viridis_d(option = "viridis")+
      labs(title = "d", x = "1 - Specificity", y = "Sensitivity", color = "Module Scoring Method") +
      theme_minimal()
      
  }
  
  # Return summary results and plot
  list_to_return = list(auc_results = res_table, roc_curves = plot, roc_curves_obj=roc_curves)
  return(list_to_return)
}



#' compute Average Silhouette Width score 
#' @param data module matrix
#' @param clustes vector with clusters labels
calculate_asw <- function(mat, clusters) {
  if (length(unique(clusters)) < 2) return(NA)
  silhouette_result <- silhouette(clusters, dist(mat))
  return(mean(silhouette_result[, 3]))
}


#' compute calinski_harabasz index (CHI)
#' @param data module matrix
#' @param clustes vector with clusters labels
# calculate_calinski_harabasz <- function(data, clusters) {
#   if (length(unique(clusters)) < 2) return(NA)
#   return(as.numeric(cluster.stats(dist(data), clusters)$ch))
# }

#' Compute Modularity of each class
#' @param mat module matrix
#' @param clustes vector with clusters label' 
#' @param k number of neighbors
calculate_modularity = function(mat, clusters, k=10){
  knn <- get.knn(mat, k = k)
  
  edges <- cbind(
    rep(1:nrow(mat), each = k),
    as.vector(knn$nn.index)
  )
  
  g <- graph_from_edgelist(edges, directed = FALSE)
  g <- simplify(g)
  
  q <- modularity(g, membership = clusters)
  return(q)
}


#' Compute Davies Bouldin index
#' @param mat module matrix
#' @param clustes vector with clusters label' 
calculate_davies_bouldin = function(mat, clusters){
  db=intCriteria(traj = as.matrix(mat),
                 part = as.integer(clusters),
                 crit = "Davies_Bouldin")
  return(db$davies_bouldin)
}

#' Compute DBCV index
#' @param mat module matrix
#' @param clustes vector with clusters label' 
calculate_dbcv_index=function(mat, clusters){
  return(dbcv_index(data=mat, partition=clusters, noiseLabel = -1))
}

#' compute ASW and CHI for one module matrix
#' @param mat_name name of the method tested 
#' @param mat module matrix for this method
evaluate_clustering <- function(mat_name, mat) {
  results_local <- data.frame()
  
  for (k in 2:min(5, nrow(mat))) {  # Limiter k à max nombre de lignes
    message(k, " clusters")
    clusters <- tryCatch({
      cutree(hclust(dist(mat)), k)
    }, error = function(e) return(rep(NA, nrow(mat))))
    
    if (any(is.na(clusters))) next  # Ignorer si clustering échoue
    
    # ch_index <- calculate_calinski_harabasz(mat, clusters)
    asw_score <- calculate_asw(mat, clusters)
    # modularity = calculate_modularity(mat, clusters)
    dbcv_index = calculate_dbcv_index(mat, clusters)
    davies_bouldin = calculate_davies_bouldin(mat, clusters)
    
    # results_local <- rbind(results_local, data.frame(
    #   Dataset = mat_name,
    #   Clusters = k,
    #   Metric = "Calinski-Harabasz",
    #   Value = ch_index
    # ))
    
    
    results_local <- rbind(results_local, data.frame(
      Dataset = mat_name,
      Clusters = k,
      Metric = "ASW Score",
      Value = asw_score
    ))
    
    results_local <- rbind(results_local, data.frame(
      Dataset = mat_name,
      Clusters = k,
      Metric = "Davies Bouldin",
      Value = davies_bouldin
    ))
    
    
    # results_local <- rbind(results_local, data.frame(
    #   Dataset = mat_name,
    #   Clusters = k,
    #   Metric = "Modularity",
    #   Value = modularity
    # ))
    
    
    results_local <- rbind(results_local, data.frame(
      Dataset = mat_name,
      Clusters = k,
      Metric = "DBCV index",
      Value = dbcv_index
    ))
    
  }
  return(results_local)
}



#' compute clustering quality evaluation on matrix list
#' @param matrix_list list with objects from gedo() or compute_module_matrix(). the list must be named.
#' @param num_cores Number of cores for parallelisation. By default, num_cores <- detectCores() - 1. 
compute_clustering_quality = function(matrix_list, num_cores=NULL){
  if(is.null(num_cores)){num_cores <- detectCores() - 1}
  
  
  
  plan(multisession, workers = num_cores)
  #plan(sequential)
  #handlers("txtprogressbar")
  
  #1. compute clustering metrics
  
  export = list(matrix_list=matrix_list,
                evaluate_clustering=evaluate_clustering,
                calculate_asw=calculate_asw
                # ,calculate_calinski_harabasz=calculate_calinski_harabasz
                )
  
  
  #Version parallélisée
  # cl_res_list=future_lapply(names(matrix_list), FUN = function(mat_name){
  #   cat(paste0(mat_name,"\n"))
  #   return(evaluate_clustering(mat_name, mat = matrix_list[[mat_name]]$module_matrix))
  # }, future.seed=T, future.globals = export, future.packages=packages_list)
  
  cl_res_list=lapply(names(matrix_list), FUN = function(mat_name){
    message(mat_name)
    return(evaluate_clustering(mat_name, mat = matrix_list[[mat_name]]$module_matrix))
  })
  
  results <- data.table(do.call(rbind, cl_res_list))
  
  
  #2. ploting results
  means <- results[, .(mean_value = mean(Value, na.rm = TRUE)), by = .(Dataset, Metric)]
  desired_order <- names(matrix_list)
  results[, Dataset := factor(Dataset, levels = desired_order)]
  means[, Dataset := factor(Dataset, levels = desired_order)]
  
  p <- ggplot(results, aes(x = Clusters, y = Value, color = Dataset)) +
    geom_line() +
    geom_point() +
    scale_colour_viridis_d(option = "viridis")+
    geom_hline(
      data = means,
      aes(yintercept = mean_value, color = Dataset),
      linetype = "dashed",
      show.legend = FALSE  
    ) +
    facet_wrap(~Metric, scales = "free_y") +
    labs(x = "Number of clusters", y = "Indice", color = "Module Scoring Method") +
    theme_minimal()
  
  
  return(list(plot=p, res=results))
  
}


#' plot PHATE projections of module matrices
#' @param data one module matrix to project : data.table with modules in column and individuals in line.
#' @param k number of neighbors for PHATE
#' @param meta_data a data.table with clinical data to merge
plot_phate=function(data, k, meta_data, omic_id){
  clinical_vars = c("diag","EXPRESSION_PRECISESADS_IFN","AUTOANTIBODY_SSA","AUTOANTIBODY_SSA_52","AUTOANTIBODY_SSA_60","AUTOANTIBODY_SSB")
  
  #1. phate
  
  set.seed(123)
  phate = data.table(phateR::phate(data = data,ndim = 2, knn.dist.method = "euclidean", mds.dist.method = "euclidean", knn=k, verbose=T, mds.solver = "smacof")$embedding)
  colnames(phate)=c("PHATE1","PHATE2")
  phate$SAMPLING_OMIC_NUMBER= omic_id
  phate[meta_data, on = "SAMPLING_OMIC_NUMBER", (clinical_vars) := mget(paste0("i.", clinical_vars))]
  phate$diag=diag
  
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
  
  #2. plots
  p1 = ggplot(phate, aes(x = PHATE1, y = PHATE2, color = factor(diag))) +
    geom_point(size=1.5, alpha=0.5) +
    theme_minimal()+
    labs(color="DIAGNOSIS")+
    scale_color_manual(values = c("Control"="green","SjD"="blue"))+
    # coord_fixed()+
    theme(legend.position = "none")+ theme_cadre
  
  
  p2 = ggplot(phate, aes(x = PHATE1, y = PHATE2, color = EXPRESSION_PRECISESADS_IFN)) +
    geom_point(size = 1.5, alpha = 0.7) +
    scale_color_gradientn(colors = c("blue", "cyan", "yellow", "red")) +  
    theme_minimal()+
    labs(color="IFN score")+
    # coord_fixed()+
    theme(legend.position = "none") + theme_cadre
  
    
  
  
  # p3 = as_grob(ggplot(phate, aes(x = PHATE1, y = PHATE2, color = AUTOANTIBODY_SSA)) +
  #   geom_point(size = 1.5, alpha = 0.7) +
  #   scale_color_gradientn(colors = c("blue", "cyan", "yellow", "red")) +  
  #   theme_minimal()+
  #   labs(color="SSA Autoantibodies"))
  # 
  # p4 = as_grob(ggplot(phate, aes(x = PHATE1, y = PHATE2, color = AUTOANTIBODY_SSA_52)) +
  #   geom_point(size = 1.5, alpha = 0.7) +
  #   scale_color_gradientn(colors = c("blue", "cyan", "yellow", "red")) +  
  #   theme_minimal()+
  #   labs(color="SSA-52 Autoantibodies"))
  # 
  # p5 = as_grob(ggplot(phate, aes(x = PHATE1, y = PHATE2, color = AUTOANTIBODY_SSA_60)) +
  #   geom_point(size = 1.5, alpha = 0.7) +
  #   scale_color_gradientn(colors = c("blue", "cyan", "yellow", "red")) +  
  #   theme_minimal()+
  #   labs(color="SSA-60 Autoantibodies"))
  # 
  # p6 = as_grob(ggplot(phate, aes(x = PHATE1, y = PHATE2, color = AUTOANTIBODY_SSB)) +
  #   geom_point(size = 1.5, alpha = 0.7) +
  #   scale_color_gradientn(colors = c("blue", "cyan", "yellow", "red")) +  
  #   theme_minimal()+
  #   labs(color="SSB Autoantibodies"))
  
  # a <- wrap_elements(full = p1)
  # b <- wrap_elements(full = p2)
  # # c <- wrap_elements(full = p3)
  # # d <- wrap_elements(full = p4)
  # # e <- wrap_elements(full = p5)
  # # f <- wrap_elements(full = p6)
  # 
  # # comb <- a/b/c/d/e/f  +
  # #   plot_layout(ncol = 3, guides = "collect") +
  # #   plot_annotation(tag_levels = "a") &
  # #   theme(
  # #     plot.tag          = element_text(size = 20, face = "bold"),
  # #     plot.tag.position = c(0, 1)
  # #   )
  # 
  # comb <- a/b  +
  #   plot_layout(ncol = 2, guides = "collect") +
  #   # plot_annotation(tag_levels = "a") &
  #   theme(
  #     plot.tag          = element_text(size = 20, face = "bold"),
  #     plot.tag.position = c(0, 1)
  #   )
  
  # #3. return combined plot
  # plot_list=list(p1,p2,p3,p4,p5,p6)
  # combined_plots <- Reduce(`+`, plot_list)
  return(list(plot=list(p1, p2), phate=phate))
}





#' compute comparison of module AUC to predict Control/Diseased, according to each config
#' @param folder_for_mgedo_obj folder where are stored .rds mgedo_objects.
compute_auc_modules_grid_search = function(folder_for_data){
  
  files = list.files(folder_for_data, full.names = T)
  
  #1. Computing AUC for each module and each element of matrix_list
  auc_results <- list()
  
  for(file in files){
    print(file)
    file_name = sub("\\.rds$", "", basename(file))

    
    obj = readRDS(file = file)
    dt = obj$module_matrix
    modules = colnames(dt)
    dt$diag = obj$diag
    
    for (pathway in setdiff(x = colnames(dt), y = c("diag"))) {
      # print(paste("method=", file_name, ". module=", pathway))
      auc_value <- auc(dt$diag, dt[[pathway]], quiet=T)
      
      auc_results <- append(auc_results, list(data.table(
        config = file_name, 
        pathway = pathway,
        auc = as.numeric(auc_value)
      )))
    }
    
  }
  dt_auc_results <- rbindlist(auc_results, use.names = TRUE, fill = TRUE)
  
  dt_auc_results[, median_auc := median(auc), by=config]
  dt_auc_results <- dt_auc_results[order(-median_auc)]
  
  
  levels=rev(unique(dt_auc_results$config))
  dt_auc_results[, config := factor(config, levels = levels)]
  

  
  #3. AUC face to face
  auc_wide <- dcast(dt_auc_results, pathway ~ config, value.var = "auc")
  
  cat("Finished !\n")


  return(auc_wide)
}






#' plot UMAP/PHATE projections of module matrices
#' @param data one module matrix to project : data.table with modules in column and individuals in line.
#' @param k number of neighbors for PHATE
#' @param meta_data a data.table with clinical data to merge
plot_module_matrix=function(data, k, meta_data, omic_id, method=c("umap","phate"), metric=c("cosine","euclidean","correlation"), decay=NULL){
  clinical_vars = c("diag","EXPRESSION_PRECISESADS_IFN","AUTOANTIBODY_SSA_CLASS","AUTOANTIBODY_SSA_52_CLASS","AUTOANTIBODY_SSA_60_CLASS","AUTOANTIBODY_SSB_CLASS")
  
  #1. UMAP
  if(method=="umap"){
  set.seed(123)
  data_r = data.table(uwot::umap(X = data,n_neighbors = k,n_components = 2, metric = metric))
  }
  
  if(method=="phate"){
    use_python("/usr/bin/python3", required = TRUE)
    set.seed(123)
    data_r = data.table(phateR::phate(data = data,ndim = 2, knn.dist.method = metric, mds.dist.method = metric, knn=k, verbose=F)$embedding, decay=decay)
    if("decay" %in% colnames(data_r)){data_r[,decay:=NULL]}
  }
  
  #PHATE
  colnames(data_r)=c("CP_1","CP_2")
  
  data_r$SAMPLING_OMIC_NUMBER= omic_id
  data_r[meta_data, on = "SAMPLING_OMIC_NUMBER", (clinical_vars) := mget(paste0("i.", clinical_vars))]
  data_r$diag=diag
  
  
  #2. plots
  p1 = as_grob(ggplot(data_r, aes(x = CP_1, y = CP_2, color = factor(diag))) +
                 geom_point(size=1.5, alpha=0.5) +
                 theme_minimal()+
                 labs(color="DIAGNOSIS")+
                 scale_color_manual(values = c("Control"="green","SjD"="blue")))+
    theme(axis.text.x=element_blank(), axis.text.y=element_blank())
  
  
  p2 = as_grob(ggplot(data_r, aes(x = CP_1, y = CP_2, color = EXPRESSION_PRECISESADS_IFN)) +
                 geom_point(size = 1.5, alpha = 0.7) +
                 scale_color_gradientn(colors = c("blue", "cyan", "yellow", "red")) +  
                 theme_minimal()+
                 labs(color="IFN score")+
                 theme(axis.text.x=element_blank(), axis.text.y=element_blank())
               
  )
  
  p3 = as_grob(ggplot(data_r, aes(x = CP_1, y = CP_2, color = AUTOANTIBODY_SSA_CLASS)) +
                 geom_point(size = 1.5, alpha = 0.7) +
                 scale_color_gradientn(colors = c("blue", "cyan", "yellow", "red")) +  
                 theme_minimal()+
                 labs(color="SSA Autoantibodies"))
  
  p4 = as_grob(ggplot(data_r, aes(x = CP_1, y = CP_2, color = AUTOANTIBODY_SSA_52_CLASS)) +
                 geom_point(size = 1.5, alpha = 0.7) +
                 scale_color_gradientn(colors = c("blue", "cyan", "yellow", "red")) +  
                 theme_minimal()+
                 labs(color="SSA-52 Autoantibodies"))
  
  p5 = as_grob(ggplot(data_r, aes(x = CP_1, y = CP_2, color = AUTOANTIBODY_SSA_60_CLASS)) +
                 geom_point(size = 1.5, alpha = 0.7) +
                 scale_color_gradientn(colors = c("blue", "cyan", "yellow", "red")) +  
                 theme_minimal()+
                 labs(color="SSA-60 Autoantibodies"))
  
  p6 = as_grob(ggplot(data_r, aes(x = CP_1, y = CP_2, color = AUTOANTIBODY_SSB_CLASS)) +
                 geom_point(size = 1.5, alpha = 0.7) +
                 scale_color_gradientn(colors = c("blue", "cyan", "yellow", "red")) +  
                 theme_minimal()+
                 labs(color="SSB Autoantibodies"))
  
  a <- wrap_elements(full = p1)
  b <- wrap_elements(full = p2)
  c <- wrap_elements(full = p3)
  d <- wrap_elements(full = p4)
  e <- wrap_elements(full = p5)
  f <- wrap_elements(full = p6)
  
  comb <- a/b/c/d/e/f  +
    plot_layout(ncol = 3, guides = "collect") +
    plot_annotation(tag_levels = "a") &
    theme(
      plot.tag          = element_text(size = 20, face = "bold"),
      plot.tag.position = c(0, 1)
    )
  
  # #3. return combined plot
  # plot_list=list(p1,p2,p3,p4,p5,p6)
  # combined_plots <- Reduce(`+`, plot_list)
  return(list(plot=comb, data_r=data_r))
}






#' plot UMAP projections of gene modules in the patient space
#' @param matrix_list list with objects from gedo() or compute_module_matrix(). the list must be named.
#' @param k the number of neighbors for UMAP. 
#' @param distance_metric the distance metric for UMAP

plot_umap_modules=function(matrix_list, k, distance_metric){
  
  #ifn modules : 
  reactome_modules <- data.table(msigdbr(species = "Homo sapiens", collection = "C7", subcollection = "IMMUNESIGDB"))
  regex_ifna <- "(?i)(interferon[-_ ]*alpha|\\bifna\\b|\\bifn[-_ ]*a\\b)"
  
  reactome_modules[,ifn_modules := grepl(regex_ifna, gs_description, perl = T)]
  ifn = unique(reactome_modules[ifn_modules==T]$gs_name)
  
  
  
  
  umap_list = list()
  for(m in names(matrix_list)){
    module_matrix = matrix_list[[m]]$module_matrix
    t_module_matrix=t(module_matrix)
    set.seed(123)
    umap = data.table(uwot::umap(X = t_module_matrix, n_neighbors = k,metric = distance_metric, n_components = 2))
    colnames(umap)=c("UMAP_1","UMAP_2")
    umap$module = rownames(t_module_matrix)
    umap[, ifn_module := module %in% ifn]
    umap_list=list.append(umap_list, umap)
  }
  names(umap_list)=names(matrix_list)
  
  
  
  plot_gedo = ggplot(umap_list$GEDO, aes(x=UMAP_1, y=UMAP_2, color=ifn_module))+
    geom_point(alpha=0.5)+ labs(title = "a")+theme_minimal()+scale_colour_manual(values = c("FALSE"="black","TRUE"="red"))
  plot_umap_gedo = ggplot(umap_list$UMAP_GEDO, aes(x=UMAP_1, y=UMAP_2, color=ifn_module))+
    geom_point(alpha=0.5)+ labs(title = "b")+theme_minimal()+scale_colour_manual(values = c("FALSE"="black","TRUE"="red"))
  plot_umap_pca1 = ggplot(umap_list$PCA1, aes(x=UMAP_1, y=UMAP_2, color=ifn_module))+
    geom_point(alpha=0.5)+ labs(title = "c")+theme_minimal()+scale_colour_manual(values = c("FALSE"="black","TRUE"="red"))
  plot_umap_m_zscore = ggplot(umap_list$MZS, aes(x=UMAP_1, y=UMAP_2, color=ifn_module))+
    geom_point(alpha=0.5)+ labs(title = "d")+theme_minimal()+scale_colour_manual(values = c("FALSE"="black","TRUE"="red"))
  
  plot_list=list(plot_gedo, plot_umap_gedo, plot_umap_pca1, plot_umap_m_zscore)
  
  combined_plots <- Reduce(`+`, plot_list)
  return(combined_plots)
  
}





#' plot UMAP projections of gene modules in the patient space
#' @param model "rf" or "knn"
#' @param k the number of neighbors for KNN. 
#' @param k_folds the number of folds 
#' @param dt_list the list of module matrices
#' @param num_cores the number of cores for parallelization
#' @param y vector with continuous values to predict

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
    
    #Load & protect data 
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
    
    #Train / test split (train only used for CV)
    set.seed(423)
    train_idx <- caret::createDataPartition(dt$y, p = 0.7, list = FALSE)
    trainData <- dt[train_idx, ]
    
    #K-fold CV on training set
    set.seed(423)
    folds <- caret::createFolds(trainData$y, k = k_folds, list = TRUE)
    
    fold_results <- lapply(seq_along(folds), function(i) {
      
      test_idx  <- folds[[i]]
      testFold <- trainData[test_idx, ]
      trainFold <- trainData[-test_idx, ]
      
      #Train model
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
      
      # Metrics 
      rmse <- sqrt(mean((obs - preds)^2))
      mae  <- mean(abs(obs - preds))
      
      list(
        fold = i,
        rmse = rmse,
        mae  = mae
      )
    })
    
    # Fold-level table 
    folds_dt <- rbindlist(lapply(fold_results, function(res) {
      data.table(
        method = dt_name,
        fold   = res$fold,
        RMSE   = res$rmse,
        MAE    = res$mae
      )
    }))
    
    # Summary table
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
  
  # Aggregate across methods
  summary_table <- rbindlist(lapply(res_list, `[[`, "summary"))
  folds_table   <- rbindlist(lapply(res_list, `[[`, "folds"))
  
  setorder(summary_table, RMSE_mean)
  
  message("Regression CV completed")
  
  return(list(
    summary = summary_table,
    folds   = folds_table
  ))
}










