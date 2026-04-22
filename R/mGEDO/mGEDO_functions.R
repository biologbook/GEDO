
# 1. General functions ----------------------------------------------------
#' unload all packages before loading usefull ones
detach_all_packages <- function() {
  basic_packages <- c("package:stats", "package:graphics", "package:grDevices",
                      "package:utils", "package:datasets", "package:methods", "package:base")
  
  loaded_packages <- search()[grepl("^package:", search())]
  to_detach <- setdiff(loaded_packages, basic_packages)
  
  for (pkg in to_detach) {
    detach(pkg, character.only = TRUE, unload = TRUE)
  }
}


#' Select CpG data for data for one gene module
#' @param Y methylomic data (CpG in column)
#' @param annotation annotation object
#' @param reactome_modules data.table with gene sets and genes
#' @param module_name nale of the gene set to extract from data (Y)

select_Y_module <- function(Y, annotation, reactome_modules, module_name) {
  #Load CpG annotation 
  data("IlluminaHumanMethylation450kanno.ilmn12.hg19")
  annotation = data.table(as.data.frame(getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)))
  
  genes_of_interest <- unique(reactome_modules[gs_name == module_name]$gene_symbol)
  
  regions_of_interest <- c("TSS1500", "TSS200", "5'UTR", "1stExon")
  
  annotation_filtered <- annotation[
    sapply(strsplit(as.character(UCSC_RefGene_Group), ";"),
           function(groups) any(groups %in% regions_of_interest)) &
      sapply(strsplit(as.character(UCSC_RefGene_Name), ";"),
             function(gene_list) any(gene_list %in% genes_of_interest))
  ]
  
  CpGs_to_keep <- annotation_filtered$Name
  
  CpGs_to_keep_in_Y <- intersect(colnames(Y), CpGs_to_keep)
  Y_filtered <- Y[, ..CpGs_to_keep_in_Y]
  
  return(Y_filtered)
} 

#' extract and save X (RNAseq) subsets for each gene set
#' @param X RNAseq data (genes in column)
#' @param modules names of the modules to treat
#' @param reactome_modules data.table with gene sets and genes
#' @param folder_for_data folder to save splited data
#' @param omic_id_x patients IDS for RNAseq
split_X=function(X=NULL,modules, reactome_modules, folder_for_data, omic_id_x=NULL){
  X_modules_dir <- "X_module_subsets_dir"
  
  if (!dir.exists(paste0(folder_for_data, X_modules_dir))) {
    dir.create(paste0(folder_for_data, X_modules_dir))
    print(paste0("Directory created to save X subsets: ", normalizePath(paste0(folder_for_data, X_modules_dir))
    ))
  }
  
  existing_files <- list.files(paste0(folder_for_data, X_modules_dir), pattern = "\\.rds$", full.names = FALSE)
  existing_modules <- gsub("^X_(.+)\\.rds$", "\\1", existing_files)
  
  modules_to_create <- setdiff(modules, existing_modules)
  
  if (length(modules_to_create) > 0) {
    library(progressr)
    handlers("txtprogressbar") 
    
    with_progress({
      p <- progressor(along = modules_to_create)
      
      invisible(lapply(modules_to_create, function(module_name) {
        p(message = sprintf("Module: %s", module_name))
        
        genes_in_module =  unique(reactome_modules[gs_name == module_name]$ensembl_gene)
        genes_in_data = colnames(X)
        genes_to_keep = genes_in_data[genes_in_data %in% genes_in_module]
        # if(length(genes_to_keep)<2){stop()}
        X_sub = X[,..genes_to_keep, with=F]
        rownames(X_sub) = omic_id_x
        
        saveRDS(X_sub, file = file.path(paste0(folder_for_data, X_modules_dir), paste0("X_", module_name, ".rds")))
      }))
    })
  } else {
    message(paste0("All the files X_*.rds are already in the folder_for_data", file.path(paste0(folder_for_data, X_modules_dir))))
  }
  return(normalizePath(paste0(folder_for_data, X_modules_dir)))
}

#' extract and save Y (CpG) subsets for each gene set
#' @param Y CpG data (CpG in column)
#' @param modules names of the modules to treat
#' @param reactome_modules data.table with gene sets and genes
#' @param folder_for_data folder to save splited data
#' @param omic_id_x patients IDS for CpG data
split_Y=function(Y=NULL,modules, reactome_modules, folder_for_data, omic_id_y=NULL){
  Y_modules_dir <- "Y_module_subsets_dir"
  
  if (!dir.exists(paste0(folder_for_data, Y_modules_dir))) {
    dir.create(paste0(folder_for_data, Y_modules_dir))
    print(paste0("Directory created to save Y subsets: ", normalizePath(paste0(folder_for_data, Y_modules_dir))
    ))
  }
  
  existing_files <- list.files(paste0(folder_for_data, Y_modules_dir), pattern = "\\.rds$", full.names = FALSE)
  existing_modules <- gsub("^Y_(.+)\\.rds$", "\\1", existing_files)
  
  modules_to_create <- setdiff(modules, existing_modules)
  
  if (length(modules_to_create) > 0) {
    library(progressr)
    handlers("txtprogressbar") 
    
    with_progress({
      p <- progressor(along = modules_to_create)
      
      invisible(lapply(modules_to_create, function(module_name) {
        p(message = sprintf("Module: %s", module_name))
        
        Y_sub <- select_Y_module(
          Y = Y,
          annotation = annotation,
          reactome_modules = reactome_modules,
          module_name = module_name
        )
        
        if(ncol(Y_sub)==0){
          cat(paste("No CpG for this gene set. Deleting ",  module_name, "in the X_ folder.\n"))
          
          file_path_module_y = file.path(paste0(folder_for_data, Y_modules_dir), paste0("Y_", module_name, ".rds"))
          file_path_module_x = gsub(x = file_path_module_y, pattern = "Y", replacement = "X")
          rm(file_path_module_x)
          
        }else{
          cat(paste(module_name, "\n"))
          cat(paste("nrow(Y_sub): ", nrow(Y_sub),"\n"))
          rownames(Y_sub) = omic_id_y
          saveRDS(Y_sub, file = file.path(paste0(folder_for_data, Y_modules_dir), paste0("Y_", module_name, ".rds")))
        }
      }))
    })
  } else {
    message(paste0("All the files Y_*.rds are already in the folder_for_data", file.path(paste0(folder_for_data, Y_modules_dir))))
  }
  return(normalizePath(paste0(folder_for_data, Y_modules_dir)))
}


#' align patient with complete subset for RNAseq and CpG, by omic_id (SAMPLING_OMIC_NUMBER)
#' @param data_module_X RNAseq subset
#' @param data_module_Y CpG subset
#' @param diag vector with diag (SjD or Control)
align_data = function(data_module_X, data_module_Y, diag){
  
  omic_id_x = rownames(data_module_X)
  omic_id_y = rownames(data_module_Y)
  
  data_module_X$id = omic_id_x
  data_module_Y$id = omic_id_y
  
  
  ids_communs <- intersect(omic_id_x, omic_id_y)
  
  diag = diag[names(diag) %in% ids_communs] #TODO : vérifier que c'est bien les bons ids ici.
  
  data_module_X <- as.data.table(data_module_X)[, id := omic_id_x]
  data_module_Y <- as.data.table(data_module_Y)[, id := omic_id_y]
  
  data_module_X <- data_module_X[id %in% ids_communs]
  data_module_Y <- data_module_Y[id %in% ids_communs]
  
  data_module_X <- data_module_X[order(match(id, ids_communs))] #ordre = ids_commmun
  data_module_Y <- data_module_Y[order(match(id, ids_communs))]
  
  data_module_X[, id := NULL]
  data_module_Y[, id := NULL]
  
  rownames(data_module_X) = ids_communs
  rownames(data_module_Y) = ids_communs
  
  
  return(list(data_module_X=data_module_X, data_module_Y=data_module_Y, diag=diag))
  
}

#' secure if function
#' @param cond condition tested
IF <- function(cond) {
  lab <- deparse(substitute(cond))
  if (!is.logical(cond) || length(cond) != 1 || is.na(cond)) {
    stop(sprintf("Bad if condition: %s (type=%s, len=%d)", 
                 lab, typeof(cond), length(cond)))
  }
  cond
}


#' Fisher test with simulated p-values
#' @param data data to test
#' @param data viable to test in data
#' @param by classes to compare
fisher_sim_test <- function(data, variable, by, ...) {
  
  tbl <- table(data[[variable]], data[[by]])
  
  fit <- fisher.test(
    tbl,
    simulate.p.value = TRUE,
    B = 50000
  )
  
  tibble::tibble(
    statistic = NA_real_,         # Fisher n’a pas de statistique
    p.value = fit$p.value,
    conf.low = NA_real_,
    conf.high = NA_real_,
    method = "Fisher (simulated)"
  )
}


# 2. mGEDO ----------------------------------------------------------------

#' Global function to run mGEDO
#' @param X RNAseq data
#' @param Y CpG data
#' @param config config file : list of all parameters needed
mgedo=function(X=NULL,Y=NULL, config){
  if(is.null(config$num_cores)){num_cores <- detectCores() - 1}else{num_cores=config$num_cores}
  
  if(config$integration_type=="XY"){
    if(config$integration_method =="nemo" & config$distance !="euclidean"){
      stop("NEMO ne fonctionne qu'avec les distances euclidiennes. Modifier config.")
    }
  }
  
  if(length(unique(config$diag))>2){
    cat("Error : more than two classes in diag\n") 
    return(NULL)}
  if(length(unique(config$diag))<2){
    cat("Error : less than two classes in diag\n") 
    return(NULL)}
  if(any(is.na(config$diag))){
    cat("Error : NA in diag\n") 
    return(NULL)}
  
  if(config$category %in% c("C7","C2")){
    reactome_modules <- data.table(msigdbr(species = "Homo sapiens", collection = config$category, subcollection = config$subcategory))
    reactome_modules <- reactome_modules[, .SD, .SDcols = c("gs_name", "ensembl_gene","gene_symbol")]
    module_number = length(unique(reactome_modules$gs_name))
    modules = unique(reactome_modules$gs_name)
  }
  if(config$category =="chaussabel"){
    chauss_modules = load_chaussabel_modules()
    reactome_modules=copy(chauss_modules)
    modules = unique(chauss_modules$gs_name)
  }
  
  
  if(config$test_setting==T){
    cat("test setting : 3 modules \n")
    modules_to_test=modules[1:3]
    reactome_modules=reactome_modules[gs_name %in% modules_to_test]
    cat(paste(length(modules_to_test), "\n"))
  }else{
    modules_to_test = copy(modules)
  }
  
  cat("Spliting X by modules\n")
  X_modules_dir=split_X(X,modules = modules_to_test, reactome_modules, folder_for_data=config$folder_for_data, config$omic_id_x)
  
  cat("Spliting Y by modules\n")
  
  Y_modules_dir=split_Y(Y,modules=modules_to_test, reactome_modules, folder_for_data=config$folder_for_data, omic_id_y = config$omic_id_y)
  
  
  if(config$integration_type=="X+Y"){
    cat("integration_type=='X+Y'\n")
    gedo_obj = compute_mgedo_X_Y(X_modules_dir=X_modules_dir,
                                 Y_modules_dir=Y_modules_dir,
                                 config=config,
                                 modules_to_test=modules_to_test)
  }
  
  if(config$integration_type=="XY"){
    cat("integration_type=='XY'\n")
    
    if(config$integration_method=="snf"){
      #ids communs
      ids_communs <- intersect(config$omic_id_x, config$omic_id_y)
      #on refait le vecteur diag pour le racourcir aux id communs, et le mettre dans le bon sens des ids.
      diag = diag[names(diag) %in% ids_communs] #TODO : vérifier que c'est bien les bons ids ici.
      #ajout au fichier config
      config$diag_filtred_by_snf = diag
    }
    
    
    gedo_obj = compute_mgedo_XY(
      X_modules_dir=X_modules_dir,
      Y_modules_dir=Y_modules_dir,
      modules_to_test=modules_to_test,
      config=config)
  }
  
  gedo_obj_to_return = list(module_matrix = gedo_obj, diag=as.factor(diag), config=config) 
  class(gedo_obj_to_return)="mGEDO"
  
  return(gedo_obj_to_return)
}



#' run UMAP to prepare data
#' @param data input data to reduce
#' @param dim_reduc_method dimension reduction method. "none" for no dimension reduction.
#' @param ncomp number of component for UMAP
#' @param dim_reduc_dist_method distance method for dimension reduction (i.g. "correlation","euclidean)
#' @param k_graph number of neighbors for UMAP
reduct_dimension=function(data, dim_reduc_method=c("none","umap"), ncomp, dim_reduc_dist_method, k_graph){
  
  if(IF(dim_reduc_method=="umap")){
    if(IF(ncol(data)<=ncomp)){
      cat("ncol(data)<=ncomp. No dimension reduction\n")
      data_reduced=copy(data)
    }else{
      # cat("umap\n")
      set.seed(123)
      data_reduced = data.table(uwot::umap(X = data, n_neighbors = k_graph, n_components = ncomp, metric = dim_reduc_dist_method,verbose = F))
      colnames(data_reduced)=paste0("V", 1:ncomp)
      rownames(data_reduced) = rownames(data)
    }
  }
  
  if(IF(dim_reduc_method=="phate")){
    if(IF(ncol(data)<=ncomp)){
      cat("ncol(data)<=ncomp. No dimension reduction\n")
      data_reduced=copy(data)
    }else{
      # cat("phate\n")
      use_python("/shared/software/miniconda/envs/r-4.4.1/bin/python", required = TRUE)
      
      set.seed(123)
      data_reduced=data.table(phateR::phate(data = data,
                                            ndim = ncomp,
                                            knn.dist.method = dim_reduc_dist_method,
                                            mds.dist.method = dim_reduc_dist_method,
                                            knn=k_graph,
                                            mds.solver = "smacof",
                                            verbose=F,
                                            n.jobs = 1
      )$embedding)
      colnames(data_reduced)=paste0("V", 1:ncomp)
      rownames(data_reduced) = rownames(data)
    }
  }
  
  if(IF(dim_reduc_method=="none")){
    data_reduced=copy(data)
  }
  
  return(data_reduced)
} 



# 3. Integration X+Y ------------------------------------------------------

#' Compute module matrices for X and Y
#' @param X_modules_dir directory for subset data (X)
#' @param Y_modules_dir directory for subset data (Y)
#' @param config config list
#' @param modules_to_test list of all modules to compute and merge
compute_mgedo_X_Y = function(X_modules_dir,
                             Y_modules_dir,
                             config,
                             modules_to_test){
  
  
  
  # cat("compute_mgedo_X_Y\n")
  
  Sys.setenv(R_PROFILE_USER = "")
  Sys.setenv(R_ENVIRON_USER = "")
  plan(multisession, workers = num_cores)
  handlers("txtprogressbar")
  # cat("mGEDO X+Y :  Parallelization on modules","\n")
  
  export_x=list(
    config=config,
    compute_graph_g_rann=compute_graph_g_rann,
    compute_ts_from_data=compute_ts_from_data,
    get_cores_from_graph=get_cores_from_graph,
    reduct_dimension=reduct_dimension,
    compute_ts_from_graph=compute_ts_from_graph,
    X_modules_dir=X_modules_dir,
    modules_to_test = modules_to_test
  )
  
  #--------------RNASEQ
  cat("compute_ts_on_modules_X\n")
  
  
  with_progress({
    p <- progressor(along = modules_to_test)
    
    X_module_matrix_list <- future_lapply(
      modules_to_test,
      function(module_name, p, config, X_modules_dir, compute_ts_from_data) {
        p(sprintf("Module: %s", module_name))
        
        X_subset <- readRDS(file.path(X_modules_dir, paste0("X_", module_name, ".rds")))
        omic_id_x = rownames(X_subset)
        #scale X_subset 
        # X_subset <- X_subset[, lapply(.SD, function(x) (x - min(x)) / (max(x) - min(x)))]
        # rownames(X_subset) = omic_id_x
        diag_after_selection = config$diag[names(config$diag) %in% omic_id_x]
        config$diag_after_selection = diag_after_selection
        if (ncol(X_subset) < 3) stop("Less than 3 genes for the module")
        
        packages_list = c("data.table", "FNN", "magrittr", "tryCatchLog","rlist",
                          "igraph", "pbapply","rgl","rdist","msigdbr","msigdbdf",
                          "future.apply", "progressr", "dplyr", "tidyr","parallel","phateR","reticulate",
                          "devtools","RANN", "BiocParallel", 
                          "IlluminaHumanMethylation450kanno.ilmn12.hg19","circlize","NEMO","SNFtool","logging","futile.logger")
        
        sapply(packages_list, require, character.only = TRUE)
        
        res = compute_ts_from_data(
          data_module = X_subset,
          config=config,
          module_name=module_name)
        return(as.numeric(res$TS))
        
      },
      p = p,
      config = config,
      X_modules_dir = X_modules_dir,
      compute_ts_from_data = compute_ts_from_data,
      future.seed = TRUE,
      future.packages = packages_list
    )
  })
  
  
  names(X_module_matrix_list)=paste0("X_",modules_to_test)
  X_module_matrix_list <- Filter(Negate(is.null), X_module_matrix_list)
  X_module_matrix <- data.table(do.call(cbind, X_module_matrix_list))
  rownames(X_module_matrix) = config$omic_id_x
  
  
  cat("TS on X computed ! \n")
  
  
  #----------METYLOMIC
  cat("compute_ts_on_modules_Y\n")
  
  
  export_y=list(
    config=config,
    compute_graph_g_rann=compute_graph_g_rann,
    compute_ts_from_data=compute_ts_from_data,
    get_cores_from_graph=get_cores_from_graph,
    reduct_dimension=reduct_dimension,
    compute_ts_from_graph=compute_ts_from_graph,
    Y_modules_dir=Y_modules_dir,
    modules_to_test = modules_to_test
  )
  
  
  with_progress({
    p <- progressor(along = modules_to_test)
    
    Y_module_matrix_list <- future_lapply(
      modules_to_test,
      function(module_name, p, config, Y_modules_dir, compute_ts_from_data) {
        p(sprintf("Module: %s", module_name))
        
        Y_subset <- readRDS(file.path(Y_modules_dir, paste0("Y_", module_name, ".rds")))
        omic_id_y = rownames(Y_subset)
        
        # Y_subset <- Y_subset[, lapply(.SD, function(x) (x - min(x)) / (max(x) - min(x)))]
        # rownames(Y_subset) = omic_id_y
        # diag_after_selection = config$diag[names(config$diag) %in% omic_id_y]
        diag_after_selection <- config$diag[omic_id_y]
        
        config$diag_after_selection = diag_after_selection
        if (ncol(Y_subset) < 3) stop("Less than 3 genes for the module")
        
        #Charge packages in worker
        packages_list = c("data.table", "FNN", "magrittr", "tryCatchLog","rlist",
                          "igraph", "pbapply","rgl","rdist","msigdbr","msigdbdf",
                          "future.apply", "progressr", "dplyr", "tidyr","parallel","phateR","reticulate",
                          "devtools","RANN", "BiocParallel", 
                          "IlluminaHumanMethylation450kanno.ilmn12.hg19","circlize","NEMO","SNFtool","logging","futile.logger")
        
        sapply(packages_list, require, character.only = TRUE)
        
        
        res = compute_ts_from_data(
          data_module = Y_subset,
          config=config,
          module_name=module_name)
        as.numeric(res$TS)
        
        
      },
      p = p,
      config = config,
      Y_modules_dir = Y_modules_dir,
      compute_ts_from_data = compute_ts_from_data,
      future.seed = TRUE,
      future.packages = packages_list
    )
  })
  
  
  
  names(Y_module_matrix_list)=paste0("Y_",modules_to_test)
  Y_module_matrix_list <- Filter(Negate(is.null), Y_module_matrix_list)
  Y_module_matrix <- data.table(do.call(cbind, Y_module_matrix_list))
  rownames(Y_module_matrix) = config$omic_id_y
  
  cat("TS on Y computed\n")
  
  
  cat("returning mgedo object\n")
  
  list_to_return=list(X_module_matrix=X_module_matrix, Y_module_matrix=Y_module_matrix)
  
  # module_matrix=cbind(X_module_matrix, Y_module_matrix)
  
  return(list_to_return)
}



#' Load BloodGen3 gene set repertory
load_chaussabel_modules = function(){
  
  url <- "https://github.com/Drinchai/BloodGen3Module/raw/master/R/sysdata.rda"
  destfile <- tempfile(fileext = ".rda")
  download.file(url, destfile, mode = "wb")
  load(destfile)
  
  ml3 <- as.data.table(Module_listGen3)
  
  setnames(ml3, old = "Gene", new = "gene_symbol")
  
  ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  
  mapping <- getBM(
    attributes = c("hgnc_symbol", "ensembl_gene_id"),
    filters = "hgnc_symbol",
    values = unique(ml3$gene_symbol),
    mart = ensembl
  )
  
  setDT(mapping)
  setnames(mapping, old = c("hgnc_symbol", "ensembl_gene_id"),
           new = c("gene_symbol", "ensembl_gene"))
  ml3_map <- merge(ml3, mapping, by = "gene_symbol", all.x = TRUE)
  ml3_final <- ml3_map[, .(gs_name = Module,
                           gene_symbol,
                           ensembl_gene,
                           Function,
                           position
  )]
  
  ml3_final <- ml3_final[!is.na(ensembl_gene)]
  n_modules <- length(unique(ml3_final$gs_name))
  
  
  
  gene_in_modules <- ml3_final[, .N, by = gene_symbol][order(-N)]
  print(gene_in_modules)
  
  total_genes <- uniqueN(ml3_final$gene_symbol)
  shared_genes <- nrow(gene_in_modules[N > 1])
  prop_shared <- shared_genes / total_genes * 100
  
  cat("Nombre total de gènes :", total_genes, "\n")
  cat("Gènes présents dans plusieurs modules :", shared_genes,
      sprintf("(%.1f%%)", prop_shared), "\n")
  
  
  hist(gene_in_modules$N,
       breaks = 50,
       main = "Distribution du nombre de modules par gène",
       xlab = "Nombre de modules distincts",
       ylab = "Nombre de gènes")
  
  return(ml3_final)
}







# 4. Integration XY -------------------------------------------------------

#' Compute module matrices for XY integration
#' @param X_modules_dir directory for subset data (X)
#' @param Y_modules_dir directory for subset data (Y)
#' @param config config list
#' @param modules_to_test list of all modules to compute and merge
compute_mgedo_XY = function(
    X_modules_dir,
    Y_modules_dir,
    modules_to_test,
    config){
  
  cat("compute_mgedo_XY\n")
  Sys.setenv(R_PROFILE_USER = "")
  Sys.setenv(R_ENVIRON_USER = "")
  plan(multisession, workers = num_cores)
  # to_test
  # plan(sequential)
  options(progressr.enable = TRUE)  
  options(future.stdout.capture = TRUE)  
  
  
  
  
  handlers("txtprogressbar")
  cat("mGEDO XY :  Parallelization on modules","\n")
  
  export=list(
    modules_to_test=modules_to_test,
    config=config,
    #functions : 
    align_data=align_data, 
    get_snf_graph=get_snf_graph, 
    get_nemo_graph=get_nemo_graph,
    reduct_dimension=reduct_dimension,
    compute_ts_XY_on_module=compute_ts_XY_on_module,
    get_cores_from_graph=get_cores_from_graph,
    compute_ts_from_graph=compute_ts_from_graph,
    X_modules_dir=X_modules_dir,
    Y_modules_dir=Y_modules_dir,
    IF=IF)
  
  
  
  cat("compute_ts_on_modules_XY\n")
  
  
  with_progress({
    p <- progressor(along = modules_to_test)
    
    XY_module_matrix_list <- future_lapply(
      modules_to_test,
      function(module_name, p, config, X_modules_dir, Y_modules_dir, compute_ts_XY_on_module) {
        p(sprintf("Module: %s", module_name))
        
        X_subset <- readRDS(file.path(X_modules_dir, paste0("X_", module_name, ".rds")))
        Y_subset <- readRDS(file.path(Y_modules_dir, paste0("Y_", module_name, ".rds")))
        
        omic_id_x = rownames(X_subset)
        omic_id_y = rownames(Y_subset)
        
        
        if (ncol(X_subset) < 3) stop("Less than 3 genes for the module")
        if (ncol(Y_subset) < 3) stop("Less than 3 CpG for the module")
        
        
        #Charge packages in worker
        packages_list = c("data.table", "FNN", "magrittr", "tryCatchLog","rlist",
                          "igraph", "pbapply","rgl","rdist","msigdbr","msigdbdf",
                          "future.apply", "progressr", "dplyr", "tidyr","parallel","phateR","reticulate",
                          "devtools","RANN", "BiocParallel", 
                          "IlluminaHumanMethylation450kanno.ilmn12.hg19","circlize","NEMO","SNFtool","logging","futile.logger")
        
        sapply(packages_list, require, character.only = TRUE)
        
        tryCatch({
          res <- compute_ts_XY_on_module(X_subset = X_subset, Y_subset = Y_subset, config=config, module_name=module_name)
          as.numeric(res$TS)
        }, error = function(e) {
          cat("Error in module:", module_name, "\nMessage:", e$message, "\n")
          NA_real_
        })
      },
      p = p,
      config = config,
      X_modules_dir = X_modules_dir,
      Y_modules_dir = Y_modules_dir,
      compute_ts_XY_on_module = compute_ts_XY_on_module,
      future.seed = TRUE,
      future.packages = packages_list
    )
  })
  
  names(XY_module_matrix_list)=modules_to_test
  XY_module_matrix_list <- Filter(Negate(is.null), XY_module_matrix_list)
  XY_module_matrix <- data.table(do.call(cbind, XY_module_matrix_list))
  
  return(XY_module_matrix)
  
}






# 5. Signatures -----------------------------------------------------------

#' Compute Mean of z scores in comparison to controls
#' @param data RNAseq or CpG data subset to score
#' @param diag vector with diagnosis (Diseased and Controls)
#' @param reference_group name of the CTRL group in diag
compute_msz=function(data, diag, reference_group){ 
  
  #2. Mean and SD of controls
  data$diag = diag
  
  gr1 = data[diag==reference_group]
  gr1[,diag:=NULL]
  data[,diag:=NULL]
  
  gr1_means = colMeans(gr1)
  gr1_sds = apply(gr1, 2, sd)
  
  #3. Z scores
  calculate_z_scores <- function(sample, means, sds) {
    (sample - means) / sds
  }
  
  z_scores <- apply(data, 1, calculate_z_scores, means = gr1_means, sds = gr1_sds)
  
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



#' Select gene list from keyword or gene
#' @param keyword keywork (single or vector) corresponding to pathway to analyze
#' @param gene_to_search gene symbol to get pathway arround this gene 
#' @param name name for the pathway
#' @param min_speci minimum sensitivity threshold to select genes
#' @param min_speci minimum specificity threshold to select genes
#' @param include_immunsigdb TRUE/FASE. If TRUE, inclure IMMUNESIGDB in gene sets analyzed
get_gene_set=function(keyword=NULL, gene_to_search=NULL, name, include_immunsigdb=F, min_sensi=0.35, min_speci=0.95){
  print(paste("Treating", name))
  
  msig_c2 <- msigdbr(species = "Homo sapiens", collection = "C2")
  msig=data.table(msig_c2)
  
  if(include_immunsigdb){
    msig_isdb <- msigdbr(species = "Homo sapiens", collection = "C7", subcollection = "IMMUNESIGDB")
    msig=data.table(rbind(msig_c2, msig_isdb))
  }
  
  gene_sets=c()
  
  if(!is.null(keyword)){
    for(k in keyword){
      gene_sets_k = unique(msig[grepl(k, gs_name)]$gs_name)  
      gene_sets=c(gene_sets, gene_sets_k)
      gene_sets=unique(gene_sets)
    }
  }
  
  if(!is.null(gene_to_search)){
    gene_sets=unique(msig[gene_symbol==gene_to_search]$gs_name)
  }
  
  
  print(paste("gene sets :"))
  print(gene_sets)
  genes = unique(msig[gs_name %in% gene_sets]$gene_symbol)
  print(paste("Nbr of genes : ", length(genes)))
  
  gene_count = msig[gs_name %in% gene_sets, .N, by="gene_symbol"]
  gene_count$gene_symbol <- factor(gene_count$gene_symbol, levels = gene_count$gene_symbol[order(gene_count$N)])
  
  plot_repr=ggplot(gene_count[N>1], aes(x=gene_symbol, y=N))+
    geom_bar(stat="identity")+
    coord_flip()+
    theme_minimal()+
    labs(title=paste("Representation of genes in", name,"gene sets"), x="gene (when > 1 gene sets)", y="Nbr of gene sets involved")
  
  other_gene_sets = unique(msig[!gs_name %in% gene_sets]$gs_name)
  
  
  # --- 1. sub-ensembles ---
  sub_pos  <- msig[gs_name %in% gene_sets]         # positive gene sets  
  sub_neg  <- msig[gs_name %in% other_gene_sets]   # negative gene sets 
  
  N_pos <- uniqueN(sub_pos$gs_name)
  N_neg <- uniqueN(sub_neg$gs_name)
  
  # --- 2. count ---
  pos_count <- sub_pos[, .(VP = uniqueN(gs_name)), by = gene_symbol]
  neg_count <- sub_neg[, .(FP = uniqueN(gs_name)), by = gene_symbol]
  
  # --- 3. Merge  ---
  res <- data.table(gene_symbol = genes)
  
  res <- merge(res, pos_count, by = "gene_symbol", all.x = TRUE)
  res <- merge(res, neg_count, by = "gene_symbol", all.x = TRUE)
  
  # Remplacer NA (gène absent) : VP = 0, FP = 0
  res[is.na(VP), VP := 0]
  res[is.na(FP), FP := 0]
  
  # --- 4. sensi / speci ---
  res[, sensi := VP / N_pos]
  res[, speci := (N_neg - FP) / N_neg]
  
  
  
  print("Plots")
  
  plot_roc=ggplot(res, aes(x=sensi, y=speci))+
    geom_point()+
    theme_minimal()+
    xlim(c(0,1))+
    ylim(c(0.9,1))+
    geom_hline(yintercept = 0.95,linetype="dashed",color="red")+
    geom_vline(xintercept = 0.35,linetype="dashed",color="red")+
    
    geom_text_repel(
      data = res[sensi > min_sensi & speci > min_speci],
      aes(label = gene_symbol),
      size = 3,
      box.padding = 0.5,
      point.padding = 0.3,
      segment.color = "grey40",
      segment.size = 0.4,
      max.overlaps = Inf
    )+
    
    # geom_text(
    #   data = res[sensi > min_sensi & speci > min_speci],
    #   aes(label = gene_symbol),
    #   hjust = -0.1, vjust = -0.1, size = 3
    # )+
    labs(x="Sensibility",y="Specificity")
  
  
  
  
  library(data.table)
  
  res[, sensi_st := (sensi - min(sensi, na.rm = TRUE)) /
        (max(sensi, na.rm = TRUE) - min(sensi, na.rm = TRUE))]
  
  
  res[, speci_st := (speci - min(speci, na.rm = TRUE)) /
        (max(speci, na.rm = TRUE) - min(speci, na.rm = TRUE))]
  
  
  
  res[, score := sensi_st * speci_st]
  res[, score_st := (score - min(score, na.rm = TRUE)) /
        (max(score, na.rm = TRUE) - min(score, na.rm = TRUE))]
  
  res$gene_symbol <- factor(res$gene_symbol, levels = res$gene_symbol[order(res$score)])
  res=res[order(-score)]
  
  
  
  #Rescaling function of score_st :
  res[, new_min := 0.5 - (score/2)]
  res[,new_max := 0.5 + (score/2)]
  
  
  
  plot_score=ggplot(res[1:50], aes(x=gene_symbol, y=score))+
    geom_bar(stat="identity")+
    coord_flip()+
    theme_minimal()+
    labs(x="best 50 genes", y="Gene Importance Score (GIS)")
  
  
  
  #final gene set
  final_gene_set = res[speci > min_speci & sensi > min_sensi]
  final_gene_set=final_gene_set[msig, ensembl_gene:=i.ensembl_gene, on="gene_symbol"]
  
  cat("numbrer of gene in final gene set : ", nrow(final_gene_set),"\n")
  
  
  #---------GSEA
  gsea_like_compute_overlap=gsea_like_compute_overlap(gene_list = final_gene_set$gene_symbol)
  
  plot_list = list(plot_roc=plot_roc, plot_score=plot_score, plot_grid=wrap_elements(full = gsea_like_compute_overlap$heatmap$gtable))   
  
  
  list_to_return=list(name=name, gene_sets_scrapped=gene_sets, genes_tested=genes, plot_list=plot_list, results=res, final_gene_set=final_gene_set, gsea_like_compute_overlap=gsea_like_compute_overlap)
  return(list_to_return)
  
}

#' Compute overlap of our gene sets with published gene sets (like in GSEA)
#' @param gene_list the list of our genes
#' @param species species to get gene sets to test for overlapping
#' @param msig_category category of gene sets to test for overlapping
#' @param msig_subcategory subcategory of gene sets to test for overlapping
#' @param min_overlap minimum number of overlapping genes for testing
#' @param fdr_cutoff p-value cut-off after FDR correction
gsea_like_compute_overlap <- function(
    gene_list,
    species = "Homo sapiens",
    msig_category = "C2",
    msig_subcategory = NULL,
    min_overlap = 2,
    fdr_cutoff = 0.05
) {
  
  
  msig <- msigdbr(
    species = species,
    collection = msig_category,
    subcollection = msig_subcategory
  ) %>%
    dplyr::select(gs_name, gene_symbol) %>%
    distinct()
  
  gene_sets <- split(msig$gene_symbol, msig$gs_name)
  
  
  
  universe_genes <- unique(msig$gene_symbol)
  
  gene_list <- intersect(gene_list, universe_genes)
  N <- length(universe_genes)
  n <- length(gene_list)
  
  
  set.seed(123)
  overlap_results <- map_df(names(gene_sets), function(gs) {
    
    gs_genes <- gene_sets[[gs]]
    K <- length(gs_genes)
    k <- length(intersect(gene_list, gs_genes))
    
    if (k < min_overlap) return(NULL)
    
    pval <- phyper(
      q = k - 1,
      m = K,
      n = N - K,
      k = n,
      lower.tail = FALSE
    )
    
    tibble(
      gene_set = gs,
      K = K,
      k = k,
      k_over_K = k / K,
      p_value = pval
    )
  })
  
  
  # 4. Correction FDR
  
  overlap_results <- overlap_results %>%
    mutate(FDR_q_value = p.adjust(p_value, method = "BH")) %>%
    arrange(FDR_q_value, p_value)
  
  significant_sets <- overlap_results %>%
    filter(FDR_q_value <= fdr_cutoff)
  
  
  
  significant_sets_heatmap = data.table(copy(significant_sets))
  significant_sets_heatmap=significant_sets_heatmap[1:10]
  
  
  heatmap_matrix <- significant_sets_heatmap$gene_set %>%
    purrr::map(~ intersect(gene_sets[[.x]], gene_list)) %>%
    tibble::enframe(name = "gene_set", value = "gene") %>%
    tidyr::unnest(gene) %>%
    dplyr::mutate(value = 1) %>%
    tidyr::pivot_wider(
      names_from = gene_set,
      values_from = value,
      values_fill = 0
    ) %>%
    tibble::column_to_rownames("gene") %>%
    as.matrix()
  colnames(heatmap_matrix)=significant_sets_heatmap$gene_set
  heatmap_matrix_t <- t(heatmap_matrix)
  
  heatmap_plot <- pheatmap::pheatmap(
    heatmap_matrix_t,
    
    cluster_rows = F,
    cluster_cols = F,
    show_rownames = TRUE,
    show_colnames = TRUE,
    
    angle_col = 90,
    fontsize_row = 8,
    fontsize_col = 9,
    
    color = c("grey80", "#08306B"),  # gris / bleu foncé GSEA-like
    breaks = c(-0.5, 0.5, 1.5),
    
    legend = TRUE,
    legend_breaks = c(0, 1),
    legend_labels = c("No overlap", "Overlap"),
    
    main = ""
  )
  
  
  list_to_return=list(
    parameters = list(
      species = species,
      category = msig_category,
      subcategory = msig_subcategory,
      gene_list_size = n,
      universe_size = N
    ),
    results_table = overlap_results,
    significant_gene_sets = significant_sets,
    heatmap_matrix = heatmap_matrix,
    heatmap = heatmap_plot
  )
  
  return(list_to_return)
}



#' Extract and save X and Y subset for signatures
#' @param gs_name name of the pathway / gene set
#' @param folder_for_res folder to save data
#' @param gene_set gene set object from get_gene_set() function
extract_data_for_signature=function(gs_name=NULL, folder_for_res, gene_set=NULL){
  
  if(!is.null(gs_name)){
    module = gs_name
    cols_to_keep = c("SAMPLING_OMIC_NUMBER", msig[gs_name==module]$ensembl_gene)
    cols_present <- intersect(cols_to_keep, colnames(rna_seq_data))
    gene_symbol = unique(msig[gs_name==module]$gene_symbol)
  }
  
  if(!is.null(gene_set)){
    module = gene_set$name
    cols_to_keep = c("SAMPLING_OMIC_NUMBER", gene_set$final_gene_set$ensembl_gene)
    cols_present <- intersect(cols_to_keep, colnames(rna_seq_data))
    gene_symbol=gene_set$final_gene_set$gene_symbol
  }
  
  if(file.exists(paste0(folder_for_res, "X_",module,".rds"))==F){
    X = rna_seq_data[, ..cols_present, with=F]
    omic_id_x = X$SAMPLING_OMIC_NUMBER
    X$SAMPLING_OMIC_NUMBER=NULL
    rownames(X) = omic_id_x
    saveRDS(X, file=paste0(folder_for_res, "X_",module,".rds"))
  }
  
  if(file.exists(paste0(folder_for_res, "Y_",module,".rds"))==F){
    Y = get_cpg_from_gene_set(Y = methylomic_data, gene_set = gene_symbol)
    
    omic_id_y = Y$SAMPLING_OMIC_NUMBER
    Y$SAMPLING_OMIC_NUMBER = NULL
    
    rownames(Y) = omic_id_y
    
    rm(rna_seq_data)
    rm(methylomic_data)
    
    saveRDS(Y, file=paste0(folder_for_res, "Y_",module,".rds"))
  }
  
}


#' Load RNAseq data in the environment
load_rna_seq_data=function(){
  cat("loading rna seq data \n")
  
  PS_brutes <- readRDS("/home/clem/GEDO/article_mGEDO/data/PS_brutes.rds")
  PS_brutes[, diag := DIAGNOSIS_DISEASE_AT_ONSET]
  PS_brutes[, control := DIAGNOSIS_ARM]
  PS_brutes[control=="Control", diag:="Control"]
  PS_brutes[diag=="SjS", diag:="SjD"]
  
  rna_seq_data = readRDS("/home/clem/GEDO/article_mGEDO/data/rna_seq_corrected_verylowvar_filtred.rds")
  # rna_seq_data[, SAMPLING_OMIC_NUMBER:=paste0("N",SAMPLING_OMIC_NUMBER)]
  rna_seq_data[PS_brutes, diag := i.DIAGNOSIS_DISEASE_AT_ONSET, on="SAMPLING_OMIC_NUMBER"]
  rna_seq_data[PS_brutes, control:=i.DIAGNOSIS_ARM, on="SAMPLING_OMIC_NUMBER"]
  rna_seq_data[control=="Control", diag:="Control"]
  rna_seq_data=rna_seq_data[,diag:=factor(diag, levels = c("Control","SjS"))]
  rna_seq_data=rna_seq_data[diag %in% c("Control","SjS")]
  rna_seq_data[diag=="SjS",diag:="SjD"]
  diag_t = rna_seq_data$diag
  diag_t=factor(diag_t,levels=c("SjD","Control"))
  
  # PS_brutes[,diag:=DIAGNOSIS_DISEASE_AT_ONSET]
  # PS_brutes[DIAGNOSIS_ARM=="Control", diag:="Control"]
  omic_id_t = rna_seq_data$SAMPLING_OMIC_NUMBER
  rna_seq_data[, diag:=NULL][, control:=NULL]
  id_diag_tab_t = data.table(SAMPLING_OMIC_NUMBER=omic_id_t, diag=diag_t)
  
  return(list(
    rna_seq_data = rna_seq_data,
    id_diag_tab_t = id_diag_tab_t,
    omic_id_t = omic_id_t
  ))
  
}

#' Load CpG methylation data in the environment
load_methylomic_data=function(){
  
  PS_brutes <- readRDS("/home/clem/GEDO/article_mGEDO/data/PS_brutes.rds")
  PS_brutes[, diag := DIAGNOSIS_DISEASE_AT_ONSET]
  PS_brutes[, control := DIAGNOSIS_ARM]
  PS_brutes[control=="Control", diag:="Control"]
  PS_brutes[diag=="SjS", diag:="SjD"]
  
  cat("Loading methylomic data...\n")
  methylomic_data <- readRDS("/home/clem/GEDO/article_mGEDO/data/t_corrected_data_BPSSSCISRB_EPIC_FINAL.rds")
  omic_id_m = methylomic_data$SAMPLING_OMIC_NUMBER
  id_diag_tab_m = data.table(SAMPLING_OMIC_NUMBER=omic_id_m, diag="")
  id_diag_tab_m[PS_brutes, diag:= i.diag, on="SAMPLING_OMIC_NUMBER"]
  id_diag_tab_m=id_diag_tab_m[SAMPLING_OMIC_NUMBER %in% omic_id_t]
  common_ids = id_diag_tab_m$SAMPLING_OMIC_NUMBER
  methylomic_data=methylomic_data[SAMPLING_OMIC_NUMBER %in% omic_id_t]
  
  return(list(methylomic_data=methylomic_data, 
              common_ids=common_ids,
              omic_id_im = omic_id_m))
}

#' Select CpG data for one gene set
#' @param Y CpG methylation data
#' @param gene_set gene set object from get_gene_set() function
get_cpg_from_gene_set = function(Y, gene_set){
  data("IlluminaHumanMethylation450kanno.ilmn12.hg19")
  annotation = data.table(as.data.frame(getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)))
  
  genes_of_interest <- gene_set
  
  regions_of_interest <- c("TSS1500", "TSS200", "5'UTR", "1stExon")
  
  annotation_filtered <- annotation[
    sapply(strsplit(as.character(UCSC_RefGene_Group), ";"),
           function(groups) any(groups %in% regions_of_interest)) &
      sapply(strsplit(as.character(UCSC_RefGene_Name), ";"),
             function(gene_list) any(gene_list %in% genes_of_interest))
  ]
  
  CpGs_to_keep <- annotation_filtered$Name
  
  CpGs_to_keep_in_Y <-  c("SAMPLING_OMIC_NUMBER", intersect(colnames(Y), CpGs_to_keep))
  Y_filtered <- Y[, ..CpGs_to_keep_in_Y]
  
  return(Y_filtered)
}


#' Compute signatures with a variety of scores and figures
#' @param config config file to compute integrated transition score
#' @param pathway gene set object from get_gene_set() function
#' @param diag vector with diag (Diseased or Control)
#' @param folder_for_data folder to get saved data subsets
#' @param integration_method "nemo" or "snf"
#' @param ponderation TRUE/FALSE. If TRUE, data are scaled according to Gene Importance Score (GIS)
#' @param scale TRUE/FALSE, If TRUE, scores is scaled between 0 and 1
compute_score_one_signature = function(config, pathway, diag, folder_for_data, integration_method, ponderation=T, scale=F){
  cat(pathway$name,"\n")
  X_subset = readRDS(paste0(folder_for_data,"/X_",pathway$name,".rds"))
  Y_subset = readRDS(paste0(folder_for_data,"/Y_",pathway$name,".rds"))
  
  omic_id_x  = rownames(X_subset)
  omic_id_y  = rownames(Y_subset)
  
  
  if(scale==T){
    X_subset <- X_subset[, lapply(.SD, function(x) (x - min(x)) / (max(x) - min(x)))]
    rownames(X_subset) = omic_id_x
    Y_subset <- Y_subset[, lapply(.SD, function(x) (x - min(x)) / (max(x) - min(x)))]
    rownames(Y_subset) = omic_id_y
  }
  
  diag_x = diag[names(diag) %in% omic_id_x]
  diag_y = diag[names(diag) %in% omic_id_y]
  
  diag_x <- diag[omic_id_x]
  diag_y <- diag[omic_id_y]
  
  if(ponderation){
    message("Ponderation by sensi-speci")
    
    #------------------For X_subset
    
    lookup <- pathway$final_gene_set[, .(ensembl_gene, new_min, new_max)]
    
    X_subset_pond = copy(X_subset)
    for (g in colnames(X_subset)) {
      
      mn_new <- lookup[ensembl_gene == g, new_min]
      mx_new <- lookup[ensembl_gene == g, new_max]
      
      old_min <- X_subset[[g]] |> min(na.rm = TRUE)
      old_max <- X_subset[[g]] |> max(na.rm = TRUE)
      
      X_subset_pond[, (g) :=
                      mn_new + ( (get(g) - old_min) * (mx_new - mn_new) / (old_max - old_min) )
      ]
    }
    
    #---------For Y_subset
    data("IlluminaHumanMethylation450kanno.ilmn12.hg19")
    annotation = data.table(as.data.frame(getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)))
    genes_of_interest <- pathway$final_gene_set$gene_symbol
    regions_of_interest <- c("TSS1500", "TSS200", "5'UTR", "1stExon")
    
    annotation_filtered <- annotation[
      sapply(strsplit(as.character(UCSC_RefGene_Group), ";"),
             function(groups) any(groups %in% regions_of_interest)) &
        sapply(strsplit(as.character(UCSC_RefGene_Name), ";"),
               function(gene_list) any(gene_list %in% genes_of_interest))
    ]
    
    
    
    
    annotation_filtered[
      , gene_symbol := {
        genes <- strsplit(UCSC_RefGene_Name, ";")
        sapply(genes, function(glist) {
          
          u <- unique(glist)
          
          
          if (length(u) == 1) return(u)
          
          
          hit <- intersect(u, genes_of_interest)
          if (length(hit) > 0) return(hit[1])   
          
          
          return(u[1])   
        })
      }
    ]
    
    
    setnames(x = annotation_filtered, old="Name", new="cpg")
    
    standardisation_table = data.table(cpg = colnames(Y_subset))
    standardisation_table[annotation_filtered, gene_symbol:=i.gene_symbol, on="cpg"]
    standardisation_table[pathway$results, new_min := i.new_min, on="gene_symbol"]
    standardisation_table[pathway$results, new_max := i.new_max, on="gene_symbol"]
    
    lookup <- standardisation_table[, .(cpg, new_min, new_max)]
    Y_subset_pond = copy(Y_subset)
    
    for (c in colnames(Y_subset)) {
      
      mn_new <- lookup[cpg == c, new_min]
      mx_new <- lookup[cpg == c, new_max]
      
      old_min <- Y_subset[[c]] |> min(na.rm = TRUE)
      old_max <- Y_subset[[c]] |> max(na.rm = TRUE)
      
      Y_subset_pond[, (c) :=
                      mn_new + ( (get(c) - old_min) * (mx_new - mn_new) / (old_max - old_min) )
      ]
    }
    X_subset_to_analyze = copy(X_subset_pond)
    Y_subset_to_analyze = copy(Y_subset_pond)
  }else{
    X_subset_to_analyze = copy(X_subset)
    Y_subset_to_analyze = copy(Y_subset)
  }
  
  
  
  #----------- X + Y 
  
  #----MSZ
  #MSZ on X
  message("MZS")
  msz_x = compute_msz(data = X_subset_to_analyze,diag = diag_x, reference_group = "Control")
  names(msz_x) = omic_id_x
  #MSZ on Y
  msz_y = compute_msz(data = Y_subset_to_analyze,diag = diag_y, reference_group = "Control")
  names(msz_y) = omic_id_y
  #----TS
  #TS on X
  message("TS X")
  ts_x_dt = compute_ts_from_data(data_module=X_subset_to_analyze, config=config, module_name=pathway$name)
  ts_x = ts_x_dt$TS 
  names(ts_x) = omic_id_x
  
  #TS on Y
  message("TS Y")
  config$diag_after_selection = diag_y
  ts_y_dt = compute_ts_from_data(data_module=Y_subset_to_analyze, config=config, module_name=pathway$name)
  ts_y = ts_y_dt$TS 
  names(ts_y) = omic_id_y
  
  
  #------------XY integrated
  
  message("DIABLO")
  
  
  #Prepare data 
  aligned_data = align_data(data_module_X = X_subset_to_analyze, data_module_Y = Y_subset_to_analyze, diag = diag)
  X_subset_complete=aligned_data$data_module_X
  Y_subset_complete=aligned_data$data_module_Y
  
  #------DIABLO on XY (to compare)
  #garder la 1ère composante principale 
  #ne pas filtrer les gènes, tt prendre dans X et dans Y
  
  dataset_list = list(trans=as.matrix(X_subset_complete), methy=as.matrix(Y_subset_complete))
  
  block_splsda_model <- block.splsda(X=dataset_list, Y=aligned_data$diag, ncomp = 1, keepX = list(trans = ncol(X_subset_complete), methy = ncol(Y_subset_complete)))
  
  # 1) Composante 1 par bloc
  comp1_trans <- block_splsda_model$variates$trans[, 1]   # RNA
  comp1_methy <- block_splsda_model$variates$methy[, 1]   # CpG
  
  # (optionnel) loadings de la composante 1
  load_trans <- block_splsda_model$loadings$trans[, 1]
  load_methy <- block_splsda_model$loadings$methy[, 1]
  
  
  
  
  
  z <- function(x) as.numeric(scale(x))
  
  #pondérer par le poids de T et M dans diablo
  w_trans <- block_splsda_model$weights$comp1[1]
  w_methy <- block_splsda_model$weights$comp1[2]
  S_raw <- w_trans * z(comp1_trans) + w_methy * z(comp1_methy)
  
  # Orienter le signe pour que "haut = activation"
  # (exemple : corriger le signe pour être positif avec la moyenne d’expression du module)
  # rna_mean <- rowMeans(X_subset)              # ou un ssGSEA/singscore de X'
  
  #orientation par la moyenne des z-scores sur x : 
  # if (cor(S_raw, msz_x, use = "complete.obs") < 0) S_raw <- -S_raw
  
  # Score final (z-scored)
  S <- z(S_raw)
  
  test_to_return_diablo = data.table(diag=diag_x, S=S)
  if(mean(test_to_return_diablo[diag=="Control"]$S) > mean(test_to_return_diablo[diag=="SjD"]$S)){
    message("returning diablo")
    S=-S
  }
  
  names(S) = rownames(X_subset_complete)
  
  #-------TS on XY
  
  message("TS XY")
  
  ts_xy_dt=compute_ts_XY_on_module(X_subset = X_subset_to_analyze,Y_subset =Y_subset_to_analyze, config=config, module_name = pathway$name) 
  ts_xy = ts_xy_dt$TS
  names(ts_xy)=omic_id_x
  
  # ts_xy_vec = ts_xy$TS
  # # rna_mean <- rowMeans(X_subset)              # ou un ssGSEA/singscore de X'
  # #orientation par msz_x
  # if (cor(ts_xy_vec, msz_x, use = "complete.obs") < 0) ts_xy_vec <- -ts_xy_vec
  
  # test_to_return_ts_xy= data.table(diag=diag_x, TS=ts_xy$TS)
  # if(mean(test_to_return_ts_xy[diag=="Control"]$TS) > mean(test_to_return_ts_xy[diag=="SjD"]$TS)){
  #   message("returning ts_xy")
  #   ts_xy[,TS := -TS]
  # }
  # 
  
  
  graph=get_nemo_graph(data_module_X=X_subset_to_analyze, data_module_Y=Y_subset_to_analyze,
                       k_graph=15,
                       k_nemo = NA)
  
  
  
  
  
  
  #rCCA and PCA
  
  # if(ponderation == T){
  #   X_subset = copy(X_subset_to_analyze)
  #   Y_subset = copy(Y_subset_to_analyze)
  #   
  #   message("PCA COMB & RCCA PONDERATED")
  # }
  
  
  message("PCA COMBINED")
  
  score_pca_combined = score_pca_combined(X_subset = X_subset_to_analyze,Y_subset=Y_subset_to_analyze)
  
  test_to_return_pca_comb = data.table(diag=diag_x, score_pca_combined=score_pca_combined)
  if(mean(test_to_return_pca_comb[diag=="Control"]$score_pca_combined) > mean(test_to_return_pca_comb[diag=="SjD"]$score_pca_combined)){
    message("returning pca combined")
    score_pca_combined=-score_pca_combined
  }
  names(score_pca_combined) = omic_id_x
  
  message("RCCA")
  
  rcca_list=score_rcca(X_subset = X_subset_to_analyze, Y_subset = Y_subset_to_analyze)
  score_rcca = rcca_list$scores_final
  test_to_return_rcca = data.table(diag=diag_x, score_rcca=score_rcca)
  if(mean(test_to_return_rcca[diag=="Control"]$score_rcca) > mean(test_to_return_rcca[diag=="SjD"]$score_rcca)){
    message("returning rcca")
    score_rcca=-score_rcca
  }
  names(score_rcca) = omic_id_x
  
  
  # scores = data.table(
  #   msz_x = msz_x,
  #   msz_y = msz_y,
  #   ts_x = ts_x,
  #   ts_y = ts_y,
  #   diablo = S,
  #   ts_xy=ts_xy,
  #   pca_comb = score_pca_combined,
  #   rcca=score_rcca
  # )
  
  library(data.table)
  
  # tous les individus présents dans au moins un vecteur
  ids <- unique(unlist(lapply(list(msz_x, msz_y, ts_x, ts_y, S, ts_xy, score_pca_combined, score_rcca), names)))
  
  scores <- data.table(
    id = ids,
    msz_x = msz_x[ids],
    msz_y = msz_y[ids],
    ts_x = ts_x[ids],
    ts_y = ts_y[ids],
    diablo = S[ids],
    ts_xy = ts_xy[ids],
    pca_comb = score_pca_combined[ids],
    rcca = score_rcca[ids]
  )
  
  #-----------Plots
  #1. PHATE des patients colorés par i. diag, ii. et + : par tous les scores calculés précédement. 
  #pb : je peux pas faire un PHATE de T+M sur les 2 datasets différents. je peux ploter le graph par contre. 
  
  
  #Correspondance ensembl - gene_symbol in X_subset
  corresp = data.table(ensembl_gene = colnames(X_subset))
  corresp[pathway$final_gene_set, gene_symbol:=i.gene_symbol, on="ensembl_gene"]
  colnames(X_subset) = as.character(corresp$gene_symbol)
  
  
  message("Plot")
  list_1 = plot_figure_graph_scores(snf_graph = graph,scores = scores, diag=diag, X_subset=X_subset)
  plot_list_1 = list_1$plot_list
  plot_list_2 = plot_rgcca_patchwork(rgcca_res = rcca_list$rcca_res, diag=rcca_list$diag_aligned)
  
  plot_list = c(plot_list_1, plot_list_2)
  
  
  message("Plots ok !")
  
  
  
  
  
  #2. combinaisons pour analyse (scatter plots)
  # TS X ~ MSZ X pour savoir si le TS va dans le sens de l'expression 
  # TS Y ~ MSZ Y pareil pour la méthylation
  
  # DIABLO ~ MSZ X pour savoir si le signe de la 1ère composante va dans le bon sens (expression)
  # DIABLO ~ MSZ Y pareil côté méthylomique
  # TS XY ~ MSZ X pareil
  # TS XY ~ MSZ Y pareil
  
  
  # list_to_return=list(
  #   scores = scores,
  #   plot = plot
  # )
  
  list_to_return=list(
    scores = scores,
    # plot_list=plot_list,
    graph = graph,
    X_subset=X_subset,
    rgcca_res = rcca_list$rcca_res,
    diag_aligned=rcca_list$diag_aligned
  )
  
  return(list_to_return)
}

#' compute figures on graph for signature detail
#' @param snf_graph snf / nemo graph 
#' @param scores data.table with scores
#' @param diag vector with diag (Diseased or Control)
#' @param X_subset RNAseq data subset
plot_figure_graph_scores <- function(snf_graph, scores, diag, X_subset) {
  
  library(data.table)
  
  mu  <- mean(scores$msz_x, na.rm = TRUE)
  sig <- sd(scores$msz_x, na.rm = TRUE)
  
  scores[, msz_x_no_out := ifelse(abs(msz_x - mu) > 3 * sig, NA, msz_x)]
  
  rng <- range(scores$msz_x_no_out, na.rm = TRUE)
  
  scores[, msz_x := (msz_x_no_out - rng[1]) / (rng[2] - rng[1])]
  
  
  set.seed(123)
  lay <- igraph::layout_with_fr(snf_graph, niter = 1200)
  
  #----------plot diag 
  
  cols_diag <- as.character(diag)
  cols_diag[cols_diag == "Control"] <- "green"
  cols_diag[cols_diag == "SjD"]     <- "blue"
  if (is.null(E(snf_graph)$weight)) E(snf_graph)$weight <- 1
  
  set.seed(123)
  p_graph <- ggplotify::as.ggplot(function() {
    op <- par(no.readonly = TRUE); on.exit(par(op), add = TRUE)
    par(mar = c(1, 1, 1, 1), xpd = NA)
    igraph::plot.igraph(
      snf_graph,
      layout = lay,
      vertex.size  = 3,
      vertex.label = NA,
      vertex.color = cols_diag,
      edge.width   = scales::rescale(E(snf_graph)$weight, to = c(0.3, 2.5)),
      main         = ""
    )
  })
  
  leg_df <- data.frame(
    lab = factor(c("Control","SjD"), levels = c("Control","SjD")),
    x = 1, y = 1
  )
  
  p_legend <- ggplot(leg_df, aes(x, y, color = lab)) +
    geom_point() +
    scale_color_manual(
      values = c(Control = "green", SjD = "blue"),
      name = "DIAGNOSIS"
    ) +
    theme_void() +
    theme(
      legend.position = "right",
      legend.title = element_text(size = 10),
      legend.text  = element_text(size = 9),
      plot.margin  = ggplot2::margin(5, 5, 5, 5, unit = "pt")
    )
  
  plot_diag <- p_graph + p_legend + patchwork::plot_layout(widths = c(1, 0.25))
  plot_diag = cowplot::as_grob(plot_diag)
  
  
  
  msz_x <- scores$msz_x; msz_y <- scores$msz_y
  ts_x  <- scores$ts_x;  ts_y  <- scores$ts_y
  diablo <- scores$diablo; ts_xy <- scores$ts_xy
  pca_comb <- scores$pca_comb; rcca <- scores$rcca
  
  safe_rng <- function(v) {
    r <- range(v, na.rm = TRUE)
    if (!is.finite(diff(r)) || diff(r) == 0) r <- r + c(-1e-6, 1e-6)
    r
  }
  rng_msz_x <- safe_rng(msz_x); rng_msz_y <- safe_rng(msz_y)
  rng_ts_x  <- safe_rng(ts_x);  rng_ts_y  <- safe_rng(ts_y)
  rng_diablo <- safe_rng(diablo); rng_ts_xy <- safe_rng(ts_xy)
  rng_pca_comb <- safe_rng(pca_comb); rng_rcca <- safe_rng(rcca)
  
  
  map_cols <- function(vals, rng) {
    cf <- scales::col_numeric(viridisLite::viridis(256),
                              domain = rng, na.color = "grey80")
    
    cf(setNames(vals, scores$SAMPLING_OMIC_NUMBER)[igraph::V(snf_graph)$name])
  }
  cols_msz_x <- map_cols(msz_x, rng_msz_x)
  cols_msz_y <- map_cols(msz_y, rng_msz_y)
  cols_ts_x  <- map_cols(ts_x,  rng_ts_x)
  cols_ts_y  <- map_cols(ts_y,  rng_ts_y)
  cols_diablo <- map_cols(diablo, rng_diablo)
  cols_ts_xy <- map_cols(ts_xy, rng_ts_xy)
  cols_pca_comb <- map_cols(pca_comb, rng_pca_comb)
  cols_rcca <- map_cols(rcca, rng_rcca)
  
  
  
  
  
  if (is.null(igraph::E(snf_graph)$weight)) igraph::E(snf_graph)$weight <- 1
  
  set.seed(123)
  gp_with_legend <- function(g, cols, title, lay, rng, viridis_option = "D") {
    # Graphe (sans légende)
    p_graph <- ggplotify::as.ggplot(function() {
      op <- par(no.readonly = TRUE); on.exit(par(op), add = TRUE)
      par(mar = c(1, 1, 1, 1), xpd = NA)
      set.seed(123)
      igraph::plot.igraph(
        g,
        layout       = lay,
        vertex.size  = 3,
        vertex.label = NA,
        vertex.color = cols,
        edge.width   = scales::rescale(igraph::E(g)$weight, to = c(0.3, 2.5)),
        main         = title
      )
    })
    
    d <- data.frame(x = 1, y = 1, val = seq(rng[1], rng[2], length.out = 10))
    p_leg <- ggplot2::ggplot(d, ggplot2::aes(x, y, color = val)) +
      ggplot2::geom_point(alpha = 0) +  # juste pour activer l’échelle
      ggplot2::scale_color_viridis_c(limits = rng, option = viridis_option, name = title) +
      ggplot2::guides(color = ggplot2::guide_colorbar(barheight = grid::unit(40, "mm"))) +
      ggplot2::theme_void() +
      ggplot2::theme(
        legend.position = "right",
        legend.title    = ggplot2::element_text(size = 9),
        legend.text     = ggplot2::element_text(size = 8),
        plot.margin     = ggplot2::margin(5, 5, 5, 5, unit = "pt")
      )
    
    plot=p_graph + p_leg + patchwork::plot_layout(widths = c(1, 0.18))
    plot=as_grob(plot)
    return(plot)
  }
  
  
  
  p_msz_x  <- gp_with_legend(snf_graph, cols_msz_x,  "",  lay, rng_msz_x)
  p_msz_y  <- gp_with_legend(snf_graph, cols_msz_y,  "",  lay, rng_msz_y)
  p_ts_x   <- gp_with_legend(snf_graph, cols_ts_x,   "",   lay, rng_ts_x)
  p_ts_y   <- gp_with_legend(snf_graph, cols_ts_y,   "",   lay, rng_ts_y)
  p_diablo <- gp_with_legend(snf_graph, cols_diablo, "", lay, rng_diablo)
  p_ts_xy  <- gp_with_legend(snf_graph, cols_ts_xy,  "",  lay, rng_ts_xy)
  
  p_pca_comb <- gp_with_legend(snf_graph, cols_pca_comb,  "pca_comb",  lay, rng_pca_comb)
  p_rcca  <- gp_with_legend(snf_graph, cols_rcca,  "rcca",  lay, rng_rcca)
  
  
  #----------Heatmap 
  
  #scale X_subset for heatmap
  omic_id_x  = rownames(X_subset)
  X_subset <- X_subset[, lapply(.SD, function(x) (x - min(x)) / (max(x) - min(x)))]
  rownames(X_subset) = omic_id_x
  
  
  clean_and_scale_matrix <- function(X_subset) {
    
    library(data.table)
    
    # --- 1. Sauvegarde des rownames
    rn <- rownames(X_subset)
    
    # Conversion en data.table (sans toucher aux rownames)
    DT <- as.data.table(X_subset)
    
    # --- 2. Détection des outliers (par colonne)
    mu  <- DT[, lapply(.SD, mean, na.rm = TRUE)]
    sig <- DT[, lapply(.SD, sd,   na.rm = TRUE)]
    
    outlier_mat <- DT[, Map(function(x, m, s) {
      if (!is.finite(s) || s == 0) {
        rep(FALSE, length(x))
      } else {
        abs(x - m) > 3 * s
      }
    }, .SD, mu, sig)]
    
    outlier_dt <- as.data.table(outlier_mat)
    
    # Lignes avec ≥ 1 outlier
    rows_with_outlier <- outlier_dt[, Reduce(`|`, .SD)]
    
    # IDs supprimés
    removed_ids <- rn[rows_with_outlier]
    
    # Filtrage
    DT_clean <- DT[!rows_with_outlier]
    rn_clean <- rn[!rows_with_outlier]
    
    # --- 3. Rescaling (0-1 par colonne)
    DT_scaled <- DT_clean[, lapply(.SD, function(x) {
      r <- range(x, na.rm = TRUE)
      if (!is.finite(diff(r)) || diff(r) == 0) {
        return(rep(0, .N))
      }
      (x - r[1]) / (r[2] - r[1])
    })]
    
    # --- 4. Réattribution des rownames
    DT_scaled <- as.data.frame(DT_scaled)
    rownames(DT_scaled) <- rn_clean
    
    # --- Output structuré
    return(list(
      X_scaled = DT_scaled,
      removed_ids = removed_ids
    ))
  }
  
  res <- clean_and_scale_matrix(X_subset)
  
  X_subset_clean_scaled <- res$X_scaled
  
  module_matrix = as.matrix(X_subset_clean_scaled)
  id=1:nrow(module_matrix)
  module_matrix_t=t(module_matrix)
  colnames(module_matrix_t)=rownames(X_subset_clean_scaled)
  
  
  
  #2. Legends
  annotation_col=data.frame(diag=diag[names(diag) %in% rownames(X_subset_clean_scaled)])
  rownames(annotation_col)=names(diag[names(diag) %in% rownames(X_subset_clean_scaled)])
  
  reference_group="Control"
  another_group = "SjD"
  
  overlap_colors <- list(
    diag = setNames(c("green", "blue"), c(reference_group, another_group))
  )
  
  
  #3. Heatmap
  cat("heatmap  \n")
  
  heatmap=pheatmap::pheatmap(
    module_matrix_t,
    annotation_col = annotation_col,   
    annotation_colors = overlap_colors,
    cluster_rows = T,               
    cluster_cols = T,               
    show_rownames = T,              
    show_colnames = FALSE,             
    scale = "none",                    
    annotation_names_row=F,
    main=""
  )
  
  diag_aligned <- diag[rownames(module_matrix)]
  
  ord <- c(
    rownames(module_matrix)[diag_aligned == "SjD"],
    rownames(module_matrix)[diag_aligned == "Control"]
  )
  
  module_matrix_ord <- module_matrix[ord, ]
  annotation_col_ord <- data.frame(
    diag = diag[rownames(X_subset_clean_scaled)]
  )
  
  rownames(annotation_col_ord) <- rownames(X_subset_clean_scaled)
  module_matrix_ord_t=t(module_matrix_ord)
  
  heatmap_ctrl_sjd=pheatmap::pheatmap(
    module_matrix_ord_t,
    annotation_col = annotation_col_ord,
    annotation_colors = overlap_colors,
    cluster_rows = TRUE,        
    cluster_cols = FALSE,       
    show_rownames = TRUE,
    show_colnames = FALSE,
    scale = "none",
    annotation_names_row = FALSE,
    main = ""
  )
  
  
  cat("heatmap ok \n")
  
  
  heatmap_grob = as_grob(heatmap$gtable)
  heatmap_sjd_ctrl_grob = as_grob(heatmap_ctrl_sjd$gtable)
  
  
  #------------Plot pairs 
  scores$diag=diag
  
  
  plot_list=list(plot_diag=plot_diag,
                 p_msz_x=p_msz_x,
                 p_msz_y=p_msz_y,
                 p_ts_x=p_ts_x, p_ts_y=p_ts_y, p_diablo=p_diablo, p_ts_xy=p_ts_xy,  p_pca_comb=p_pca_comb, p_rcca=p_rcca, heatmap_grob=heatmap_grob, heatmap_sjd_ctrl_grob=heatmap_sjd_ctrl_grob)
  
  
  return(list(scores=scores, plot_list=plot_list))
}

#' compute rCCA figures for signature detail
#' @param rcca_res rCCA model trained (two components)
#' @param comp_plot UMAP/PHATE axis to test for correlation in circle plot
#' @param circle_cutoff r² threshold for correlation circle plot
#' @param legend TRUE/FASE. Legend for plots
#' @param diag vector with diag (Diseased or Control)
plot_rgcca_patchwork <- function(rgcca_res,
                                 # comp_bar = 1,
                                 comp_plot = c(1,2),
                                 circle_cutoff = 0.5,
                                 #imgCor_cutoff = 0.7,
                                 legend = TRUE,
                                 diag) {
  
  message("===== START plot_rgcca_patchwork =====")
  flush.console()
  
  require(mixOmics)
  require(ggplot2)
  require(patchwork)
  library(grid)
  library(gridGraphics)
  
  
  ################################
  ### 1. Circle plot
  ################################
  
  message("[1/3] Building circle plot...")
  
  p1 <- wrap_elements(
    full = grid.grabExpr({
      plotVar(
        rgcca_res,
        comp = comp_plot,
        cutoff = circle_cutoff,
        legend = legend,
        title = "Variable correlation circle",
        style = "ggplot2",
        plot = TRUE
      )
    })
  )
  message("✔ p1 OK")
  
  
  ################################
  ### 2. imgCor combined
  ################################
  
  message("[2/3] Building imgCor combined...")
  
  # 
  #   tmp <- tempfile(fileext = ".png")
  #   
  #   png(tmp, width = 2000, height = 2000, res = 300)
  #   
  #   imgCor(
  #     X = rgcca_res$X$RNA,
  #     Y = rgcca_res$X$CpG,
  #     type = "combine",
  #     X.var.names = TRUE,
  #     Y.var.names = TRUE,
  #     title = TRUE,
  #     symkey = TRUE
  #   )
  #   
  #   dev.off()
  #   
  #   p2 <- wrap_elements(full = grid::rasterGrob(png::readPNG(tmp)))
  
  
  
  
  cor_mat <- cor(
    rgcca_res$X$RNA,
    rgcca_res$X$CpG,
    use = "pairwise.complete.obs"
  )
  
  
  library(pheatmap)
  
  #lim <- max(abs(cor_mat), na.rm = TRUE)
  #breaks <- seq(-lim, lim, length.out = 101)
  breaks <- seq(-1, 1, length.out = 101)
  
  
  p2 <- wrap_elements(
    full = grid.grabExpr({
      
      pheatmap(
        cor_mat,
        cluster_rows = TRUE,
        cluster_cols = TRUE,
        scale = "none",
        color = colorRampPalette(c("blue", "white", "red"))(100),
        breaks=breaks,
        border_color = NA,
        main = "",
        fontsize = 5
      )
      
    })
  )
  message("✔ p2 OK")
  
  ################################
  ### 3. Arrow plot
  ################################
  
  message("[3/3] Building arrow plot...")
  
  # p3 <- wrap_elements(
  #   full = grid.grabExpr({
  #     plotIndiv(
  #       rgcca_res,
  #       comp = comp_plot,
  #       arrow = TRUE,
  #       legend = legend,
  #       title = "Individuals projection (arrow plot)"
  #     )
  #     gridGraphics::grid.echo()
  #   })
  # )
  
  
  
  
  
  diag <- factor(diag, levels = c("SjD", "Control"))
  
  cols <- c("SjD" = "green", "Control" = "blue")
  
  p3 <- wrap_elements(
    full = grid.grabExpr({
      plotIndiv(
        rgcca_res,
        comp = comp_plot,
        group = diag,
        col = cols,
        ind.names = FALSE,   # remplace les IDs par des points
        arrow = TRUE,
        legend = TRUE,
        title = ""
      )
      gridGraphics::grid.echo()
    })
  )
  
  
  message("✔ p3 OK")
  message("===== SUCCESS plot_rgcca_patchwork =====")
  
  plot_list=list(corr_circle_plot=p1, corr_heatmap=p2, arrow_plot=p3)
  
  return(plot_list)
}



#' compute rCCA scores and rCCA two component object
#' @param X_subset X data subset
#' @param Y_subset Y data subset
#' @param ncomp ncomp in rCCA model for scoring
score_rcca <- function(X_subset, Y_subset, ncomp = 1) {
  
  library(mixOmics)
  omic_id_x = rownames(X_subset)
  
  common_ids <- intersect(rownames(X_subset), rownames(Y_subset))
  rna_only_ids <- setdiff(rownames(X_subset), common_ids)
  
  common_ids <- intersect(rownames(X_subset), rownames(Y_subset))
  
  aligned_data = align_data(data_module_X = X_subset, data_module_Y = Y_subset, diag = diag)
  X_subset_complete=aligned_data$data_module_X
  Y_subset_complete=aligned_data$data_module_Y
  
  Xc = scale(X_subset_complete)
  Yc = scale(Y_subset_complete)
  
  
  center_X <- attr(Xc, "scaled:center")
  scale_X  <- attr(Xc, "scaled:scale")
  
  center_Y <- attr(Yc, "scaled:center")
  scale_Y  <- attr(Yc, "scaled:scale")
  
  # 3. Multi-block model
  # blocks <- list(RNA = Xc, CpG = Yc)
  design <- matrix(c(0,1,
                     1,0), 2, 2, byrow = TRUE)
  
  rgcca_res <- wrapper.rgcca(
    X = list(RNA = Xc, CpG = Yc),
    design = design,
    ncomp = ncomp
  )
  
  #Traduction in gene symbol
  mapping = data.table(ensembl_gene=colnames(Xc))
  mapping[msig, gene_symbol:=i.gene_symbol, on="ensembl_gene"]
  colnames(Xc)=mapping$gene_symbol
  rgcca_res_ncomp2 <- wrapper.rgcca(
    X = list(RNA = Xc, CpG = Yc),
    design = design,
    ncomp = 2
  )
  
  t_complete <- rgcca_res$variates$RNA[,1]
  
  # 5. Loadings RNA
  loadings_rna <- rgcca_res$loadings$RNA[,1]
  
  
  X_subset_mat = as.matrix(X_subset)
  rownames(X_subset_mat)=omic_id_x
  
  
  # 6. Projection RNA-only
  if (length(rna_only_ids) > 0) {
    
    X_rna_only_scaled <- scale(
      X_subset_mat[rna_only_ids, , drop = FALSE],
      center = center_X,
      scale  = scale_X
    )
    
    t_rna_only <- as.vector(X_rna_only_scaled %*% loadings_rna)
    
  } else {
    t_rna_only <- numeric(0)
  }
  
  # 7. Combinaison
  scores_raw <- numeric(nrow(X_subset_mat))
  names(scores_raw) <- rownames(X_subset_mat)
  
  scores_raw[common_ids] <- t_complete
  scores_raw[rna_only_ids] <- t_rna_only
  
  # 8. Standardisation 
  mu <- mean(t_complete)
  sdv <- sd(t_complete)
  
  scores_final <- (scores_raw - mu) / sdv
  
  return(list(scores_final=scores_final, rcca_res=rgcca_res_ncomp2, diag_aligned=aligned_data$diag))
}


#' compute PCA combined method for unsupervised scoring
#' @param X_subset X data subset
#' @param Y_subset Y data subset
score_pca_combined <- function(X_subset, Y_subset) {
  
  # Align IDs
  ids_all <- union(rownames(X_subset), rownames(Y_subset))
  
  # --- RNA PCA ---
  X_scaled <- scale(X_subset)
  pca_X <- prcomp(X_scaled, center = FALSE, scale. = FALSE)
  
  PCX <- rep(NA, length(ids_all))
  names(PCX) <- ids_all
  PCX[rownames(X_subset)] <- pca_X$x[,1]
  PCX <- scale(PCX)
  
  varX <- pca_X$sdev[1]^2
  
  # --- CpG PCA ---
  Y_scaled <- scale(Y_subset)
  pca_Y <- prcomp(Y_scaled, center = FALSE, scale. = FALSE)
  
  PCY <- rep(NA, length(ids_all))
  names(PCY) <- ids_all
  PCY[rownames(Y_subset)] <- pca_Y$x[,1]
  PCY <- scale(PCY)
  
  varY <- pca_Y$sdev[1]^2
  
  w1 <- varX / (varX + varY)
  w2 <- varY / (varX + varY)
  
  score <- w1 * PCX + w2 * PCY
  
  # fallback if Y missing 
  missing_Y <- is.na(PCY)
  score[missing_Y] <- PCX[missing_Y]
  
  return(as.numeric(score))
}



#' Compute signature differences after rescaling
#' @param metrics metrics to analyze (ex : rcca)
#' @param prefix prefix for new variables
compute_signature_differences <- function(
    scores_cd40,
    scores_il2,
    scores_ifn,
    metrics,
    prefix = "norm_diff"
){
  
  eps <- 1e-6
  
  # =============================
  # 1. Rename columns
  # =============================
  
  rename_with_suffix <- function(dt, suffix){
    cols <- setdiff(colnames(dt), "SAMPLING_OMIC_NUMBER")
    setnames(dt, old = cols, new = paste0(cols, "_", suffix))
    return(dt)
  }
  
  scores_cd40 <- rename_with_suffix(copy(scores_cd40), "CD40")
  scores_il2  <- rename_with_suffix(copy(scores_il2),  "IL2")
  scores_ifn  <- rename_with_suffix(copy(scores_ifn),  "IFN")
  
  # =============================
  # 2. Merge all scores
  # =============================
  
  scores_all <- Reduce(function(x,y)
    merge(x,y, by="SAMPLING_OMIC_NUMBER", all=TRUE),
    list(scores_cd40, scores_il2, scores_ifn))
  
  scores_scaled <- copy(scores_all)
  
  # =============================
  # 3. Scaling per comparison
  # =============================
  
  for(m in metrics){
    
    col_cd40 <- paste0(m, "_CD40")
    col_ifn  <- paste0(m, "_IFN")
    col_il2  <- paste0(m, "_IL2")
    
    # CD40 vs IFN (ref CD40)
    mu <- mean(scores_all[[col_cd40]], na.rm=TRUE)
    sdv <- sd(scores_all[[col_cd40]], na.rm=TRUE)
    
    scores_scaled[, paste0(m, "_CD40_scaled_cd40_ifn") :=
                    (get(col_cd40) - mu)/(sdv+eps)]
    scores_scaled[, paste0(m, "_IFN_scaled_cd40_ifn") :=
                    (get(col_ifn) - mu)/(sdv+eps)]
    
    # IL2 vs IFN (ref IL2)
    mu <- mean(scores_all[[col_il2]], na.rm=TRUE)
    sdv <- sd(scores_all[[col_il2]], na.rm=TRUE)
    
    scores_scaled[, paste0(m, "_IL2_scaled_il2_ifn") :=
                    (get(col_il2) - mu)/(sdv+eps)]
    scores_scaled[, paste0(m, "_IFN_scaled_il2_ifn") :=
                    (get(col_ifn) - mu)/(sdv+eps)]
    
    # CD40 vs IL2 (ref CD40)
    mu <- mean(scores_all[[col_cd40]], na.rm=TRUE)
    sdv <- sd(scores_all[[col_cd40]], na.rm=TRUE)
    
    scores_scaled[, paste0(m, "_CD40_scaled_cd40_il2") :=
                    (get(col_cd40) - mu)/(sdv+eps)]
    scores_scaled[, paste0(m, "_IL2_scaled_cd40_il2") :=
                    (get(col_il2) - mu)/(sdv+eps)]
  }
  
  # =============================
  # 4. Simple differences
  # =============================
  
  for(m in metrics){
    
    scores_scaled[, paste0("ratio_", m, "_cd40_ifn") :=
                    get(paste0(m, "_CD40_scaled_cd40_ifn")) -
                    get(paste0(m, "_IFN_scaled_cd40_ifn"))]
    
    scores_scaled[, paste0("ratio_", m, "_il2_ifn") :=
                    get(paste0(m, "_IL2_scaled_il2_ifn")) -
                    get(paste0(m, "_IFN_scaled_il2_ifn"))]
    
    scores_scaled[, paste0("ratio_", m, "_cd40_il2") :=
                    get(paste0(m, "_CD40_scaled_cd40_il2")) -
                    get(paste0(m, "_IL2_scaled_cd40_il2"))]
  }
  
  
  
  return(scores_scaled)
}


#' Function to aggregate a list with signatures plots
#' @param x list of signatures results
#' @param diag vector with diag (Diseased or Control)
get_figure_for_signature = function(x, diag){
  list_1 = plot_figure_graph_scores(snf_graph = x$graph,
                                    scores = x$scores,
                                    diag=diag,
                                    X_subset=x$X_subset)
  
  plot_list_1 = list_1$plot_list
  
  plot_list_2 = plot_rgcca_patchwork(rgcca_res = x$rgcca_res,
                                     diag=x$diag_aligned)
  
  plot_list = c(plot_list_1, plot_list_2)
  
  return(plot_list)
}



# 6. Plots ----------------------------------------------------------------

#' Preparation of a list with mgedo objects and metadata
#' @param folder_for_mgedo_obj folder to get mGEDO data
prepare_mgedo_obj_metadata = function(folder_for_mgedo_obj){
  #Get mgedo_obj
  files = list.files(path = folder_for_mgedo_obj,full.names = T)
  list_mgedo_obj = lapply(files,readRDS)
  
  names= gsub(x = files, pattern = folder_for_mgedo_obj, replacement = "")
  names=gsub(x=names,pattern = ".rds",replacement = "")
  names(list_mgedo_obj) = names
  
  setkey(PS_brutes, SAMPLING_OMIC_NUMBER)
  
  ifn_score_communs <- PS_brutes[ids_communs, EXPRESSION_PRECISESADS_IFN]
  
  ifn_score <- PS_brutes[omic_id_t, EXPRESSION_PRECISESADS_IFN]
  
  meta_data = list(
    id_diag_tab_t=id_diag_tab_t,
    id_diag_tab_m=id_diag_tab_m,
    ids_communs=ids_communs,
    diag_t = diag_t,
    diag_m = diag_m,
    diag=diag,
    ifn_score_communs=ifn_score_communs,
    ifn_score=ifn_score
  )
  
  return(list(list_mgedo_obj=list_mgedo_obj, meta_data=meta_data))
}


#' Function to plot and save UMAP and PHATE embeddings from mGEDO objects
#' @param mgedo_obj_list_metadata result from prepare_mgedo_obj_metadata()
#' @param folder_for_data folder to get mGEDO objects
#' @param folder_for_res folder to save embeddings
#' @param umap_k k for umap
#' @param phate_k k for phate
#' @param metric distance metric ("correlation","euclidean","cosine")
#' @param min_dist min_dist for umap
#' @param spread spread for umap
#' @param repulsion_strength repulsion_strength for umap 
#' @param negative_sample_rate negative_sample_rate for umap 
#' @param n_epochs n_epochs for umap 
#' @param decay decay for phate 
#' @param t t for phate
#' @param do_pca TRUE/FALSE. If TRUE, UMAP and PHATE include PCA before embedding
#' @param pca_ncomp Number of PCA components if do_pca
plot_phate_umap_on_mgdedo=function(mgedo_obj_list_metadata,
                                   folder_for_data,
                                   folder_for_res,
                                   umap_k=10,
                                   phate_k = 12,
                                   metric="cosine",
                                   min_dist = 0.001,
                                   spread = 1,
                                   repulsion_strength =1,
                                   negative_sample_rate = 10,
                                   n_epochs=500,
                                   decay = 60,
                                   t=10,
                                   do_pca=T,
                                   pca_ncomp=40
){
  
  
  module_matrix_list = mgedo_obj_list_metadata$list_mgedo_obj
  meta_data_list = mgedo_obj_list_metadata$meta_data
  
  plot_list=list()
  embedding_list = list()
  
  for(i in names(module_matrix_list)){
    print(paste(gene_sets_to_test, ": ", i))
    module_matrix = module_matrix_list[[i]]$module_matrix
    
    if(inherits(module_matrix, "list")){
      message("X+Y integration")
      X_module_matrix = module_matrix$X_module_matrix
      Y_module_matrix = module_matrix$Y_module_matrix
      
      use_python("/shared/software/miniconda/envs/r-4.4.1/bin/python", required = TRUE)
      set.seed(123)
      X_phate_2d = data.table(phateR::phate(data = X_module_matrix,ndim = 2, knn.dist.method = metric, mds.dist.method = metric, knn=phate_k, verbose=F, decay=decay, mds.solver = "smacof",t=t,n.landmark = NULL,npca = pca_ncomp)$embedding)
      set.seed(123)
      Y_phate_2d = data.table(phateR::phate(data = Y_module_matrix,ndim = 2, knn.dist.method = metric, mds.dist.method = metric, knn=phate_k, verbose=F, decay=decay, mds.solver = "smacof", t=t,n.landmark = NULL,npca = pca_ncomp)$embedding)
      set.seed(123)
      X_umap_2d = data.table(uwot::umap(X = X_module_matrix,
                                        n_neighbors = umap_k,
                                        n_components = 2,
                                        metric = metric,
                                        min_dist=min_dist,
                                        spread = spread,
                                        repulsion_strength =repulsion_strength,
                                        negative_sample_rate =negative_sample_rate,
                                        n_epochs=n_epochs,
                                        pca = pca_ncomp,n_threads = 1))
      colnames(X_umap_2d)= c("UMAP1","UMAP2")
      set.seed(123)
      Y_umap_2d = data.table(uwot::umap(X = Y_module_matrix,
                                        n_neighbors = umap_k,
                                        n_components = 2,
                                        metric = metric,
                                        min_dist=min_dist,
                                        spread = spread,
                                        repulsion_strength =repulsion_strength,
                                        negative_sample_rate =negative_sample_rate,
                                        n_epochs=n_epochs,
                                        pca = pca_ncomp,n_threads = 1))
      colnames(Y_umap_2d)= c("UMAP1","UMAP2")
      
      
      #For X 
      X_omic_id =  meta_data_list$id_diag_tab_t$SAMPLING_OMIC_NUMBER
      X_ifn_score = meta_data_list$ifn_score
      X_diag = meta_data_list$diag
      
      X_phate_2d$SAMPLING_OMIC_NUMBER = X_omic_id
      X_phate_2d$diag = X_diag
      X_phate_2d$EXPRESSION_PRECISESADS_IFN = X_ifn_score
      
      X_umap_2d$SAMPLING_OMIC_NUMBER = X_omic_id
      X_umap_2d$diag = X_diag
      X_umap_2d$EXPRESSION_PRECISESADS_IFN = X_ifn_score
      
      
      #For Y 
      Y_umap_2d[, SAMPLING_OMIC_NUMBER := rownames(Y_module_matrix)]
      Y_umap_2d <- merge(
        Y_umap_2d,
        PS_brutes[, .(SAMPLING_OMIC_NUMBER, diag, EXPRESSION_PRECISESADS_IFN)],
        by = "SAMPLING_OMIC_NUMBER",
        all.x = TRUE
      )
      
      
      Y_phate_2d[, SAMPLING_OMIC_NUMBER := rownames(Y_module_matrix)]
      Y_phate_2d <- merge(
        Y_phate_2d,
        PS_brutes[, .(SAMPLING_OMIC_NUMBER, diag, EXPRESSION_PRECISESADS_IFN)],
        by = "SAMPLING_OMIC_NUMBER",
        all.x = TRUE
      )
      
      
      X_plot_phate = plot_phate_2d(X_phate_2d, title = paste(gene_sets_to_test,"-X-", i))
      Y_plot_phate = plot_phate_2d(Y_phate_2d, title = paste(gene_sets_to_test,"-Y-", i))
      
      X_plot_umap = plot_umap_2d(X_umap_2d, title = paste(gene_sets_to_test,"-X-", i))
      Y_plot_umap = plot_umap_2d(Y_umap_2d, title = paste(gene_sets_to_test,"-Y-", i))
      
      embedding_list =  list.append(embedding_list, X_umap_2d=X_umap_2d, Y_umap_2d=Y_umap_2d, X_phate_2d=X_phate_2d, Y_phate_2d =Y_phate_2d)
      plot_list=list.append(plot_list, X_plot_phate=X_plot_phate, Y_plot_phate=Y_plot_phate, X_plot_umap=X_plot_umap, Y_plot_umap=Y_plot_umap)
      
      
    }
    
    
    if("data.frame" %in% class(module_matrix)){
      message("XY integration")
      # print(paste("dim(module_matrix):", dim(module_matrix) ))
      
      if(anyNA(module_matrix)){
        na_cols <- names(module_matrix)[sapply(module_matrix, function(x) any(is.na(x)))]
        print(paste(gene_sets_to_test, ": ", i,". cols with NA :", na_cols))
        print(paste(gene_sets_to_test, ": ", i,". nbr of columns :", length(na_cols), "_ Deleting columns."))
        module_matrix[, (na_cols) := NULL]
      }
      
      
      print(paste(gene_sets_to_test, ": ", i," - computing phate..."))
      use_python("/shared/software/miniconda/envs/r-4.4.1/bin/python", required = TRUE)
      set.seed(123)
      phate_2d = data.table(phateR::phate(data = module_matrix,ndim = 2, knn.dist.method = metric, mds.dist.method = metric, knn=phate_k, verbose=F, decay=decay, mds.solver = "smacof",t=t,n.landmark = NULL,npca = pca_ncomp)$embedding)
      
      
      print(paste(gene_sets_to_test, ": ", i," - computing umap..."))
      set.seed(123)
      umap_2d = data.table(uwot::umap(X = module_matrix,
                                      n_neighbors = umap_k,
                                      n_components = 2,
                                      metric = metric,
                                      min_dist=min_dist,
                                      spread = spread,
                                      repulsion_strength =repulsion_strength,
                                      negative_sample_rate =negative_sample_rate,
                                      n_epochs=n_epochs,
                                      pca = pca_ncomp,n_threads = 1))
      colnames(umap_2d)= c("UMAP1","UMAP2")
      
      if(module_matrix_list[[i]]$config$integration_method == "snf"){
        # cat("snf \n")
        omic_id =  meta_data_list$ids_communs
        ifn_score = meta_data_list$ifn_score_communs
        diag = meta_data_list$diag_m
      }
      
      if(module_matrix_list[[i]]$config$integration_method == "nemo"){
        # cat("nemo \n")
        omic_id =  meta_data_list$id_diag_tab_t$SAMPLING_OMIC_NUMBER
        ifn_score = meta_data_list$ifn_score
        diag = meta_data_list$diag
      }
      
      
      phate_2d$SAMPLING_OMIC_NUMBER = omic_id
      phate_2d$diag = diag
      phate_2d$EXPRESSION_PRECISESADS_IFN = ifn_score
      
      umap_2d$SAMPLING_OMIC_NUMBER = omic_id
      umap_2d$diag = diag
      umap_2d$EXPRESSION_PRECISESADS_IFN = ifn_score
      
      embedding_list =  list.append(embedding_list,umap_2d = umap_2d, phate_2d = phate_2d)
      
      plot_phate = plot_phate_2d(phate_2d, title = paste("PHATE-", gene_sets_to_test,"-", i, "-", module_matrix_list[[i]]$config$integration_method))
      plot_umap = plot_umap_2d(umap_2d, title = paste("UMAP-", gene_sets_to_test,"-", i, "-", module_matrix_list[[i]]$config$integration_method))
      
      # ggsave(plot, filename= paste0(folder_for_res,"phate_module_matrix_", i, "_",gene_sets_to_test,".png"))
      
      
      plot_list=list.append(plot_list, plot_phate=plot_phate, plot_umap=plot_umap)
    }
    
  }
  
  return(list(embedding_list=embedding_list, plot_list = plot_list))
}




#' Create plots on PHATE embedding
#' @param phate data.table with phate axis and metadata
#' @param title title for plots
plot_phate_2d = function(phate, title=""){
  
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
  
  ac= c( "AUTOANTIBODY_SSA_52_CLASS",  "AUTOANTIBODY_SSA_60_CLASS" ,"AUTOANTIBODY_SSB_CLASS")
  
  phate[PS_brutes, AUTOANTIBODY_SSA_52_CLASS:=i.AUTOANTIBODY_SSA_52_CLASS, on="SAMPLING_OMIC_NUMBER"]
  phate[PS_brutes, AUTOANTIBODY_SSA_60_CLASS:=i.AUTOANTIBODY_SSA_60_CLASS, on="SAMPLING_OMIC_NUMBER"]
  phate[PS_brutes, AUTOANTIBODY_SSB_CLASS:=i.AUTOANTIBODY_SSB_CLASS, on="SAMPLING_OMIC_NUMBER"]
  phate[, (ac):=lapply(.SD, as.factor), .SDcols = ac]
  
  p_diag = ggplot(phate, aes(x = PHATE1, y = PHATE2, color = factor(diag))) +
    geom_point(size=1.5, alpha=0.5) +
    theme_minimal()+
    labs(color="DIAGNOSIS")+
    scale_color_manual(values = c("Control"="green","SjD"="blue"))+
    theme(axis.text.x=element_blank(), axis.text.y=element_blank())+
    theme_cadre
  
  p_kirou = ggplot(phate, aes(x = PHATE1, y = PHATE2, color = EXPRESSION_PRECISESADS_IFN)) +
    geom_point(size = 1.5, alpha = 0.7) +
    scale_color_viridis_c(option = "H") +
    theme_minimal()+
    labs(color="IFN Kirou score")+theme_cadre
  
  
  
  
  
  p_ssa52 = ggplot(phate, aes(x = PHATE1, y = PHATE2, color = AUTOANTIBODY_SSA_52_CLASS)) +
    geom_point(size = 1.5, alpha = 0.7) +
    scale_color_manual(values=c("0"="#0D0887FF", "1"="#6A00A8FF","2"="#B12A90FF", "3"="#E16462FF", "4"="#FCA636FF"))+
    theme_minimal()+
    labs(color="SSA-52 Autoantibodies")+theme(axis.text.x=element_blank(), axis.text.y=element_blank())+theme_cadre
  
  p_ssa60 =ggplot(phate, aes(x = PHATE1, y = PHATE2, color = AUTOANTIBODY_SSA_60_CLASS)) +
    geom_point(size = 1.5, alpha = 0.7) +
    scale_color_manual(values=c("0"="#0D0887FF", "1"="#6A00A8FF","2"="#B12A90FF", "3"="#E16462FF", "4"="#FCA636FF"))+
    theme_minimal()+
    labs(color="SSA-60 Autoantibodies")+theme(axis.text.x=element_blank(), axis.text.y=element_blank())+theme_cadre
  
  p_ssb = ggplot(phate, aes(x = PHATE1, y = PHATE2, color = AUTOANTIBODY_SSB_CLASS)) +
    geom_point(size = 1.5, alpha = 0.7) +
    scale_color_manual(values=c("0"="#0D0887FF", "1"="#6A00A8FF","2"="#B12A90FF", "3"="#E16462FF", "4"="#FCA636FF"))+
    theme_minimal()+
    labs(color="SSB Autoantibodies")+theme(axis.text.x=element_blank(), axis.text.y=element_blank())+theme_cadre
  
  plot_list = list(p_diag=p_diag, p_kirou=p_kirou, p_ssa52=p_ssa52, p_ssa60=p_ssa60, p_ssb=p_ssb)
  
  return(plot_list)
  
}




#' Create plots on UMAP embedding
#' @param umap data.table with umap axis and metadata
#' @param title title for plots
plot_umap_2d = function(umap, title=""){
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
  
  
  ac= c( "AUTOANTIBODY_SSA_52_CLASS",  "AUTOANTIBODY_SSA_60_CLASS" ,"AUTOANTIBODY_SSB_CLASS")
  
  umap[PS_brutes, AUTOANTIBODY_SSA_52_CLASS:=i.AUTOANTIBODY_SSA_52_CLASS, on="SAMPLING_OMIC_NUMBER"]
  umap[PS_brutes, AUTOANTIBODY_SSA_60_CLASS:=i.AUTOANTIBODY_SSA_60_CLASS, on="SAMPLING_OMIC_NUMBER"]
  umap[PS_brutes, AUTOANTIBODY_SSB_CLASS:=i.AUTOANTIBODY_SSB_CLASS, on="SAMPLING_OMIC_NUMBER"]
  umap[, (ac):=lapply(.SD, as.factor), .SDcols = ac]
  
  p_diag = ggplot(umap, aes(x = UMAP1, y = UMAP2, color = factor(diag))) +
    geom_point(size=1.5, alpha=0.5) +
    theme_minimal()+
    labs(color="DIAGNOSIS")+
    scale_color_manual(values = c("Control"="green","SjD"="blue"))+ theme(axis.text.x=element_blank(), axis.text.y=element_blank())+theme_cadre
  
  p_kirou =ggplot(umap, aes(x = UMAP1, y = UMAP2, color = EXPRESSION_PRECISESADS_IFN)) +
    geom_point(size = 1.5, alpha = 0.7) +
    scale_color_gradientn(colors = c("blue", "cyan", "yellow", "red")) +  
    theme_minimal()+
    labs(color="IFN Kirou score")+theme_cadre
  
  p_ssa52 = ggplot(umap, aes(x = UMAP1, y = UMAP2, color = AUTOANTIBODY_SSA_52_CLASS)) +
    geom_point(size = 1.5, alpha = 0.7) +
    scale_color_manual(values=c("0"="#0D0887FF", "1"="#6A00A8FF","2"="#B12A90FF", "3"="#E16462FF", "4"="#FCA636FF"))+
    theme_minimal()+
    labs(color="SSA-52 Autoantibodies")+theme(axis.text.x=element_blank(), axis.text.y=element_blank())+theme_cadre
  
  p_ssa60 = ggplot(umap, aes(x = UMAP1, y = UMAP2, color = AUTOANTIBODY_SSA_60_CLASS)) +
    geom_point(size = 1.5, alpha = 0.7) +
    scale_color_manual(values=c("0"="#0D0887FF", "1"="#6A00A8FF","2"="#B12A90FF", "3"="#E16462FF", "4"="#FCA636FF"))+
    theme_minimal()+
    labs(color="SSA-60 Autoantibodies")+theme(axis.text.x=element_blank(), axis.text.y=element_blank())+theme_cadre
  
  p_ssb = ggplot(umap, aes(x = UMAP1, y = UMAP2, color = AUTOANTIBODY_SSB_CLASS)) +
    geom_point(size = 1.5, alpha = 0.7) +
    scale_color_manual(values=c("0"="#0D0887FF", "1"="#6A00A8FF","2"="#B12A90FF", "3"="#E16462FF", "4"="#FCA636FF"))+
    theme_minimal()+
    labs(color="SSB Autoantibodies")+theme(axis.text.x=element_blank(), axis.text.y=element_blank())+theme_cadre
  
  plot_list = list(p_diag=p_diag, p_kirou=p_kirou, p_ssa52=p_ssa52, p_ssa60=p_ssa60, p_ssb=p_ssb)
  
  return(plot_list)
  
}


#' Create correlation circle plot
#' @param module_matrix module matrix 
#' @param embedding phate or umap data.table
#' @param dim_names names of umap and phate axis
#' @param method method for correlation
#' @param top_n_labels number of biological functions to detail
#' @param point_size size of points
corr_circle_plot <- function(module_matrix,
                             embedding,
                             dim_names = c("Dim1", "Dim2"),
                             method = "pearson",
                             top_n_labels = 50,
                             point_size = 3) {
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
  require(data.table)
  require(ggplot2)
  require(dplyr)
  library(ggrepel )
  
  if (ncol(embedding) != 2)
    stop("embedding must contain exactly 2 numeric columns")
  
  if (nrow(module_matrix) != nrow(embedding))
    stop("module_matrix and embedding must have same number of rows")
  
  module_matrix <- as.data.frame(module_matrix)
  embedding <- as.data.frame(embedding)
  
  module_matrix[] <- lapply(module_matrix, as.numeric)
  embedding[] <- lapply(embedding, as.numeric)
  
  # --- Correlations ---
  cor_mat <- sapply(module_matrix, function(mod) {
    c(
      cor(mod, embedding[[1]], method = method, use = "pairwise.complete.obs"),
      cor(mod, embedding[[2]], method = method, use = "pairwise.complete.obs")
    )
  })
  
  cor_df <- data.frame(
    module = colnames(module_matrix),
    x = cor_mat[1, ],
    y = cor_mat[2, ]
  )
  
  # --- Distance to origin ---
  cor_df$distance <- sqrt(cor_df$x^2 + cor_df$y^2)
  
  # --- Selection top modules ---
  cor_df <- cor_df %>%
    arrange(desc(distance))
  
  cor_df$label <- ""
  cor_df$label[1:min(top_n_labels, nrow(cor_df))] <-
    cor_df$module[1:min(top_n_labels, nrow(cor_df))]
  
  
  
  highlight_set <- c("Interferon", "Type 1 Interferon")
  
  cor_df$IFN <- ifelse(
    cor_df$module %in% highlight_set,
    "IFN or Type I IFN",
    "Other"
  )
  cor_df=data.table(cor_df)
  cor_df[, IFN:=factor(IFN, levels=c("Other","IFN or Type I IFN"))]
  
  top_n <- min(top_n_labels, nrow(cor_df))
  
  top_modules <- cor_df[1:top_n, ]
  
  label_df <- top_modules[!duplicated(top_modules$module), ]
  
  circle <- data.frame(
    x = cos(seq(0, 2*pi, length.out = 500)),
    y = sin(seq(0, 2*pi, length.out = 500))
  )
  
  circle_inner <- data.frame(
    x = 0.5 * cos(seq(0, 2*pi, length.out = 500)),
    y = 0.5 * sin(seq(0, 2*pi, length.out = 500))
  )
  
  
  
  # --- Plot ---
  p <- ggplot() +
    
    geom_path(data = circle, aes(x, y), color = "black") +
    geom_path(data = circle_inner, aes(x, y), color = "black") +
    
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    
    geom_point(
      data = cor_df[IFN=="Other"],
      aes(x = x, y = y, color = IFN),
      size = point_size, alpha=0.7
    )+
    geom_point(
      data = cor_df[IFN=="IFN or Type I IFN"],
      aes(x = x, y = y, color = IFN),
      size = point_size,
      alpha=0.7
    )+
    
    scale_color_manual(
      values = c(
        "IFN or Type I IFN" = "red",
        "Other" = "grey70"
      )
    )+
    
    geom_text_repel(
      data = label_df,
      aes(x = x, y = y, label = module),
      size = 3,
      box.padding = 0.5,
      point.padding = 0.3,
      segment.color = "grey40",
      segment.size = 0.4,
      max.overlaps = Inf
    )+
    
    # coord_equal(xlim = c(-1, 1), ylim = c(-1, 1)) +
    theme_minimal() +
    labs(
      x = dim_names[1],
      y = dim_names[2],
      title = ""
    )+theme_cadre
  
  return(p)
}


#' Plot clinical variables (bars)
#' @param dt clinical data to test
#' @param vars variables to test
#' @param class names of autoantibody variables
#' @param cluster_col groups to compare
#' @param color_col diagnosis vector
#' @param ncol number of columns for figure
plot_significant_vars <- function(dt, vars, class, cluster_col = "cluster", color_col = "diag",ncol=4) {
  
  dt <- as.data.table(dt)
  
  vars_quanti <- vars[sapply(dt[, ..vars], is.numeric)]
  vars_quali  <- vars[!vars %in% vars_quanti]
  
  # ⚠️ convertir character → factor pour barplots
  for(v in vars_quali){
    if(is.character(dt[[v]])) dt[[v]] <- factor(dt[[v]])
  }
  
  # ---- 1. BOXPLOTS  ----
  
  comparisons <- combn(levels(dt$cd40_ifn_class), 2, simplify = FALSE)
  
  plots_quanti <- lapply(vars_quanti, function(v) {
    
    df_v = dt[!is.na(dt[[v]])]
    df_v = df_v[!is.na(df_v[[cluster_col]])]
    
    if (all(df_v$diag == "SjD")) {
      df_v[, diag := factor(diag, levels = "SjD")]
      df_v[, (cluster_col) := droplevels(get(cluster_col))]
    }
    comparisons <- combn(levels(df_v$cd40_ifn_class), 2, simplify = FALSE)
    
    ggplot(df_v, aes(x = get(cluster_col), y = .data[[v]]))+
      # pas de remplissage
      geom_jitter(aes(color = diag),                            
                  width = 0.15, size = 1.6, alpha = 0.8) +
      geom_boxplot(fill = NA, color = "grey25") + 
      # scale_color_manual(values=c("SjD"="blue","CTRL"="green"))+
      scale_color_manual(values=c("SjD"="blue","CTRL"="green"))+
      
      guides(fill = "none") +
      theme_minimal()+
      labs(title=v, y="",x="CD40-IFN class", color="DIAGNOSIS")+
      stat_compare_means(
        comparisons = comparisons,
        method      = "wilcox.test",
        label       = "p.signif",
        hide.ns     = F,
        na.rm=T
      )+
      theme(
        plot.title = element_text(size = 10)  
      )+
      coord_flip()
    
  })
  
  # ---- 2. BARPLOTS ----
  plots_quali <- lapply(vars_quali, function(v) {
    
    dt2 <- dt[, .SD, .SDcols = c(cluster_col, v)][complete.cases(dt[, .SD, .SDcols = c(cluster_col, v)])]
    df_plot <- dt2[
      , .(n = .N), 
      by = c(cluster_col, v)
    ][
      , prop := n / sum(n), by = cluster_col
    ]
    
    
    df_plot[[v]] <- factor(
      df_plot[[v]],
      levels = c(setdiff(levels(df_plot[[v]]), "Unknown"), "Unknown")
    )
    df_plot[[v]] <- forcats::fct_explicit_na(df_plot[[v]], na_level = "NA")
    
    n_totals <- df_plot[, .(n_total = sum(n)), by = cd40_ifn_class]
    
    
    if(v %in% class){
      
      ggplot(df_plot, aes_string(x = cluster_col, y = "prop", fill = v)) +
        geom_col(position = "fill") +
        geom_text(
          data = n_totals,
          aes(x = cd40_ifn_class,
              y = 1.05,
              label = paste0("n=", n_total)),
          inherit.aes = FALSE,
          size = 4,
          # fontface = "bold",
          color = "black"
        ) +
        scale_y_continuous(
          labels = scales::percent,
          expand = expansion(mult = c(0, 0.1))  
        ) +
        theme_minimal(base_size = 14) +
        labs(title = v, x="CD40-IFN class", y = "Proportion", fill = NULL)+
        scale_fill_manual(values=c("0"="#0D0887FF", "1"="#6A00A8FF","2"="#B12A90FF", "3"="#E16462FF", "4"="#FCA636FF","Unknown"="grey30"))+
        theme_minimal(base_size = 14)+
        theme(
          plot.title = element_text(size = 10)  # ← taille du titre
        )+
        coord_flip(clip = "off")
      
      
    }else{
      if("Low" %in% unique(df_plot[[v]]) | "Moderate" %in% unique(df_plot[[v]]) | "High" %in% unique(df_plot[[v]])){
        
        df_plot[[v]] = factor(df_plot[[v]], levels=c("No","Low","Moderate"))
        
        ggplot(df_plot, aes_string(x = cluster_col, y = "prop", fill = v)) +
          geom_col(position = "fill") +
          geom_text(
            data = n_totals,
            aes(x = cd40_ifn_class,
                y = 1.05,
                label = paste0("n=", n_total)),
            inherit.aes = FALSE,
            size = 4,
            # fontface = "bold",
            color = "black"
          ) +
          scale_y_continuous(
            labels = scales::percent,
            expand = expansion(mult = c(0, 0.1)) 
          ) +
          theme_minimal(base_size = 14) +
          labs(title = v, x="CD40-IFN class", y = "Proportion", fill = NULL)+
          scale_fill_manual(
            values = c(
              "Low"      = "#6A00A8FF",  
              "Moderate" = "#FCA636FF",  
              "High"     = "yellow",
              "No"       = "black",
              "Unknown"  = "grey30",
              "NA"       = "grey80"
            )
          )+
          theme(
            plot.title = element_text(size = 10)  
          )+
          coord_flip()
        
        
      }else{
        if(v =="diag"){
          c("diag \n")
          ggplot(df_plot, aes_string(x = cluster_col, y = "prop", fill = v)) +
            geom_col(position = "fill") +
            geom_text(
              data = n_totals,
              aes(x = cd40_ifn_class,
                  y = 1.05,
                  label = paste0("n=", n_total)),
              inherit.aes = FALSE,
              size = 4,
              # fontface = "bold",
              color = "black"
            ) +
            scale_y_continuous(
              labels = scales::percent,
              expand = expansion(mult = c(0, 0.1))  
            ) +
            theme_minimal(base_size = 14) +
            labs(title = v,x="CD40-IFN class", y = "Proportion", fill = NULL)+
            scale_fill_manual(
              values = c(
                "SjD"      = "blue",  
                "CTRL" = "green")
            )+
            theme(
              plot.title = element_text(size = 10)  
            )+
            coord_flip()
          
          
          
        }else{
          
          levels_v <- levels(df_plot[[v]])
          
          default_pal <- scales::hue_pal()(length(levels_v))
          names(default_pal) <- levels_v
          if ("Unknown" %in% levels_v) {
            default_pal["Unknown"] <- "grey30"
          }
          
          
          ggplot(df_plot, aes_string(x = cluster_col, y = "prop", fill = v)) +
            geom_col(position = "fill") +
            geom_text(
              data = n_totals,
              aes(x = cd40_ifn_class,
                  y = 1.05,
                  label = paste0("n=", n_total)),
              inherit.aes = FALSE,
              size = 4,
              # fontface = "bold",
              color = "black"
            ) +
            scale_y_continuous(
              labels = scales::percent,
              expand = expansion(mult = c(0, 0.1))  
            ) +
            labs(title = v,x="CD40-IFN class", y = "Proportion", fill = NULL) +
            theme_minimal(base_size = 14) +
            scale_fill_manual(values = default_pal)+
            theme(
              plot.title = element_text(size = 10) 
            )+
            coord_flip()
          
        }
      }
    }
    
  })
  
  # ---- 3. PATCHWORK ----
  
  all_plots <- c(plots_quanti, plots_quali)
  
  combined <- wrap_plots(all_plots, ncol = ncol)
  
  return(combined)
}


#' Plot clinical variables 2 (barplot)
#' @param dt clinical data to test
#' @param vars variables to test
#' @param class names of autoantibody variables
#' @param cluster_col groups to compare
#' @param color_col diagnosis vector
#' @param ncol number of columns for figure
plot_significant_vars_2 <- function(dt, vars, class, cluster_col = "cluster", color_col = "diag", ncol=4) {
  
  dt <- as.data.table(dt)
  
  vars_quanti <- vars[sapply(dt[, ..vars], is.numeric)]
  vars_quali  <- vars[!vars %in% vars_quanti]
  
  
  for(v in vars_quali){
    if(is.character(dt[[v]])) dt[[v]] <- factor(dt[[v]])
  }
  
  # ---- 1. BOXPLOTS  ----
  
  comparisons <- combn(levels(dt$cd40_ifn_class), 2, simplify = FALSE)
  
  plots_quanti <- lapply(vars_quanti, function(v) {
    
    df_v <- dt[!is.na(dt[[v]]) & !is.na(dt[[cluster_col]])]
    
    ggplot(df_v, aes(x = get(cluster_col), y = .data[[v]]))+
      geom_jitter(aes(color = diag),
                  width = 0.15, size = 1.6, alpha = 0.8) +
      geom_boxplot(fill = NA, color = "grey25") + 
      scale_color_manual(values=c("SjD"="blue","CTRL"="green"))+
      guides(fill = "none") +
      theme_minimal()+
      labs(title=v, y="",x="CD40-IFN class")+
      stat_compare_means(
        comparisons = comparisons,
        method      = "wilcox.test",
        label       = "p.signif",
        hide.ns     = FALSE,
        na.rm=TRUE
      )+
      theme(plot.title = element_text(size = 10))+
      coord_flip()
    
  })
  
  # ---- 2. BARPLOTS + CI ----
  
  plots_quali <- lapply(vars_quali, function(v) {
    
    dt2 <- dt[, .SD, .SDcols = c(cluster_col, v)][
      complete.cases(dt[, .SD, .SDcols = c(cluster_col, v)])
    ]
    
    
    clusters <- sort(unique(dt[[cluster_col]]))
    levels_v <- sort(unique(dt[[v]]))
    
    
    grid <- CJ(
      cluster_tmp = clusters,
      value_tmp   = levels_v,
      unique = TRUE
    )
    setnames(grid, c("cluster_tmp","value_tmp"), c(cluster_col, v))
    
    
    df_plot <- dt2[, .(n = .N), by = c(cluster_col, v)]
    
    # join + remplissage
    df_plot <- grid[df_plot, on = c(cluster_col, v)]
    df_plot[is.na(n), n := 0]
    
    
    df_plot[, n_total := sum(n), by = cluster_col]
    df_plot[, prop := fifelse(n_total > 0, n / n_total, 0)]
    
    # IC Wilson
    df_plot <- df_plot[, {
      if(n_total[1] == 0){
        .(n, n_total, prop, lower = 0, upper = 0)
      } else {
        ci <- binom::binom.confint(n, n_total, method = "wilson")
        .(n, n_total, prop, lower = ci$lower, upper = ci$upper)
      }
    }, by = c(cluster_col, v)]
    
    
    df_plot[[cluster_col]] <- factor(df_plot[[cluster_col]], levels = clusters)
    df_plot[[v]] <- factor(df_plot[[v]], levels = levels_v)
    
    # ---- position dodge  ----
    
    n_levels <- length(levels_v)
    width <- 0.8
    
    df_plot[, idx := as.numeric(get(v))]
    df_plot[, x_num := as.numeric(get(cluster_col))]
    
    df_plot[, x_pos := x_num + (idx - (n_levels + 1)/2) * (width / n_levels)]
    
    # ---- n total  ----
    
    n_totals <- unique(df_plot[, .(cluster = get(cluster_col), n_total)])
    n_totals[, x_num := as.numeric(factor(cluster, levels = clusters))]
    
    
    # ---- plot ----
    
    p <- ggplot(df_plot, aes(x = x_pos, y = prop, fill = .data[[v]])) +
      
      geom_col(width = width / n_levels) +
      
      geom_errorbar(
        aes(ymin = lower, ymax = upper),
        width = 0.15
      ) +
      
      geom_text(
        data = n_totals,
        aes(x = x_num, y = 1.05, label = paste0("n=", n_total)),
        inherit.aes = FALSE,
        size = 4
      )
    
    
    
    if(v %in% class){
      
      p <- p + scale_fill_manual(values=c(
        "0"="grey30", "1"="#6A00A8FF","2"="#B12A90FF",
        "3"="#E16462FF", "4"="#FCA636FF","Unknown"="grey30"
      ))
      
    } else if(any(c("Low","Moderate","High") %in% levels_v)){
      
      df_plot[[v]] <- factor(df_plot[[v]], levels=c("No","Low","Moderate","High"))
      
      p <- p + scale_fill_manual(
        values =c(
          "Low"      = "#6A00A8FF",  
          "Moderate" = "#FCA636FF",  
          "High"     = "yellow",
          "No"       = "grey30",
          "Unknown"  = "grey30",
          "NA"       = "grey80"
        )
      )
      
      
    } else if(v == "diag") {
      
      p <- p + scale_fill_manual(values = c("SjD"="blue","CTRL"="green"))
      
    } else {
      
      default_pal <- scales::hue_pal()(length(levels_v))
      names(default_pal) <- levels_v
      
      if ("No" %in% levels_v) {
        default_pal["No"] <- "grey30"
      }
      
      p <- p + scale_fill_manual(values = default_pal)
    }
    
    
    
    p +
      scale_x_continuous(
        breaks = seq_along(clusters),
        labels = clusters
      ) +
      scale_y_continuous(
        labels = scales::percent,
        expand = expansion(mult = c(0, 0.15))
      ) +
      labs(title = v, x="CD40-IFN class", y = "Proportion", fill = NULL) +
      theme_minimal(base_size = 14) +
      theme(plot.title = element_text(size = 10)) +
      coord_flip(clip = "off")
    
  })
  
  # ---- 3. COMBINAISON ----
  
  all_plots <- c(plots_quanti, plots_quali)
  combined <- patchwork::wrap_plots(all_plots, ncol = ncol)
  
  return(combined)
}


#' Plot clinical variables (scatter plots)
#' @param dt clinical data to test
#' @param vars variables to test
#' @param class names of autoantibody variables
#' @param color_col diagnosis vector
#' @param ncol number of columns for figure
plot_significant_vars_point <- function(dt, vars, class, 
                                        color_col = "diag",
                                        ncol = 4) {
  
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(forcats)
  library(scales)
  
  
  # ---- split variables ----
  vars_quanti <- vars[sapply(dt[, ..vars], is.numeric)]
  vars_quali  <- vars[!vars %in% vars_quanti]
  
  for(v in vars_quali){
    if(is.character(dt[[v]])) dt[[v]] <- factor(dt[[v]])
  }
  
  
  # ---- 1. SCATTER QUANTI ----
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
  plots_quanti <- lapply(vars_quanti, function(v) {
    
    df_v <- dt[!is.na(get(v)) & !is.na(PHATE1) & !is.na(PHATE2)]
    
    ggplot(df_v, aes(x = PHATE1, y = PHATE2)) +
      geom_point(aes(color = .data[[v]]), size = 1.2, alpha = 0.6) +
      scale_color_viridis_c(option = "plasma", na.value = "grey80") +
      theme_minimal() +
      labs(title = v, color = v) +
      theme(
        plot.title = element_text(size = 10)
      )+theme_cadre
  })
  
  # ---- 2. SCATTER QUALI----
  plots_quali <- lapply(vars_quali, function(v) {
    
    df_v <- dt[!is.na(get(v)) & !is.na(PHATE1) & !is.na(PHATE2)]
    
    df_v[[v]] <- forcats::fct_explicit_na(df_v[[v]], na_level = "NA")
    if("Unknown" %in% levels(df_v[[v]])){
      df_v[[v]] <- factor(
        df_v[[v]],
        levels = c(setdiff(levels(df_v[[v]]), "Unknown"), "Unknown")
      )
    }
    
    if(v %in% class){
      
      pal <- c(
        "0"="#0D0887FF", 
        "1"="#6A00A8FF",
        "2"="#B12A90FF", 
        "3"="#E16462FF", 
        "4"="#FCA636FF",
        "Unknown"="grey30",
        "NA"="grey80"
      )
      
    } else if(v == "diag"){
      
      pal <- c(
        "SjD" = "blue",
        "CTRL" = "green",
        "Unknown"="grey30",
        "NA"="grey80"
      )
      
    } else if(any(c("Low","Moderate","High") %in% levels(df_v[[v]]))){
      
      df_v[[v]] <- factor(df_v[[v]], levels=c("No","Low","Moderate","High"))
      
      pal <- c(
        "Low"      = "#6A00A8FF",  
        "Moderate" = "#FCA636FF",  
        "High"     = "yellow",
        "No"       = "black",
        "Unknown"  = "grey30",
        "NA"       = "grey80"
      )
      
    } else {
      
      levels_v <- levels(df_v[[v]])
      pal <- scales::hue_pal()(length(levels_v))
      names(pal) <- levels_v
      
      if ("Unknown" %in% levels_v) pal["Unknown"] <- "grey30"
      if ("NA" %in% levels_v) pal["NA"] <- "grey80"
    }
    
    ggplot(df_v, aes(x = PHATE1, y = PHATE2)) +
      geom_point(aes(color = .data[[v]]), size = 1.2, alpha = 0.6) +
      scale_color_manual(values = pal, na.value = "grey80") +
      theme_minimal() +
      labs(title = v, color = v) +
      theme(
        plot.title = element_text(size = 10)
      )+theme_cadre
    
  })
  
  # ---- COMBINE ----
  all_plots <- c(plots_quanti, plots_quali)
  combined <- wrap_plots(all_plots, ncol = ncol)
  
  return(combined)
}


#' Function to plot correlation with published signatures
#' @param manual_signature signatures computed from home made gene sets
#' @param pub_signature signature computed from publioshed gene sets
plot_corr_signature = function(manual_signature, pub_signature){
  cor_test <- cor.test(scores_all_with_diff[[pub_signature]], scores_all_with_diff[[manual_signature]], method = "pearson")
  
  r  <- cor_test$estimate
  r2 <- r^2
  pval <- cor_test$p.value
  
  plot=ggplot(scores_all_with_diff, aes(x = .data[[manual_signature]], y = .data[[pub_signature]])) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "lm", se = FALSE, color = "red") +
    labs(
      title = paste0(
        " | r = ", round(r, 3),
        " | R² = ", round(r2, 3),
        " | p = ", signif(pval, 3)
      ),
      x = manual_signature,
      y = pub_signature
    )+
    # coord_fixed() +
    theme_bw() 
  # +
  #   theme(
  #     aspect.ratio = 1
  #   )
  return(plot)
}












