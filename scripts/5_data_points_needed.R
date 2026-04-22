
# 1. Packages and scripts  ------------------------------------------------------------

folder_script_to_source = "/home/clem/GEDO/R/"
source(file = paste0(folder_script_to_source,"GEDO.R"))
source(file = paste0(folder_script_to_source,"functions_article.R"))

detach_all_packages()
packages_list = c("data.table", "FNN", "magrittr", "tryCatchLog","rlist",
                  "igraph","pbapply","rgl","rdist","msigdbr","msigdbdf",
                  "ggplot2", "patchwork","pROC", "RANN","dbscan", "cowplot")
sapply(packages_list, require, character.only = TRUE)

folder_for_res = "/home/clem/GEDO/results/"


# 2. Data -----------------------------------------------------------------
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
diag=factor(diag, levels=c("Control","SjD"))
length(diag)


collection="C7"
subcollection="IMMUNESIGDB"

module_name="GSE1740_UNSTIM_VS_IFNA_STIMULATED_MCSF_DERIVED_MACROPHAGE_DN"
reactome_modules <- data.table(msigdbr(species = "Homo sapiens", collection = collection, subcollection = subcollection))
genes_in_module =  unique(reactome_modules[gs_name == module_name]$ensembl_gene)
genes_in_data = colnames(rna_seq_data)
genes_to_keep = genes_in_data[genes_in_data %in% genes_in_module]
data_module = rna_seq_data[,..genes_to_keep, with=F]
dim(data_module)
rm(rna_seq_data)



# 3. Function  ------------------------------------------------------------

#Computing TS : aller chercher une fonction profonde dans GEDO.R
compute_transition_score_from_data = function(data_module, diag, reference_group="Control", k_graph=5, distance="euclidean",scale_ts=F){
  
  res=data.table(diag=diag)
  res$id = 1:nrow(res)
  
  #3. Sep gr1 and gr2
  data_module$diag = diag
  groupe_1 = data_module[diag==reference_group]
  groupe_1[,diag:=NULL]
  groupe_2 = data_module[diag != reference_group]
  groupe_2[,diag:=NULL]
  data_module[, diag:=NULL]
  
  #4. LOF on gr1 and gr2
  # lof_gr1 = lof(groupe_1, minPts=k_lof)
  # lof_gr2 = lof(groupe_2, minPts=k_lof)
  # res[, lof := ifelse(diag == reference_group, lof_gr1, lof_gr2)]
  # res[diag == reference_group, core_gr1 := lof_gr1 <= quantile(lof_gr1, core_pct)]
  # res[diag != reference_group, core_gr2 := lof_gr2 <= quantile(lof_gr2, core_pct)]
  
  res[diag == reference_group, core_gr1 := T]
  res[diag != reference_group, core_gr2 := T]
  
  graph_g=compute_graph_g_rann(X = data_module, k = k_graph,distance = distance)
  

 
  core_gr1_ids = res[core_gr1==T]$id
  core_gr2_ids = res[core_gr2==T]$id
  d_core_gr1 = apply(igraph::distances(graph_g, v = res$id, to = core_gr1_ids, mode = "all", algorithm = "dijkstra"),1,mean)
  d_core_gr2 = apply(igraph::distances(graph_g, v = res$id, to = core_gr2_ids, mode = "all", algorithm = "dijkstra"),1,mean)
  res[,':='(mean_dist_core_gr1=d_core_gr1, mean_dist_core_gr2=d_core_gr2)]
  
 
  res[, TS := mean_dist_core_gr1 / (mean_dist_core_gr1 + mean_dist_core_gr2)]
  
 
  if(scale_ts==T){
    res[, TS := (TS - min(TS)) / (max(TS) - min(TS))]
  }
  return(res)
}
res_all=compute_transition_score_from_data(data_module = data_module, diag=diag)
res_all[,diag:=as.character(diag)]
auc(res_all$diag, res_all$TS)

library(data.table)
library(pROC)

granularity = seq(from = 10, to = 600, by = 10)

test_ts_robustness <- function(
    data_module,
    diag,
    batch_sizes = seq(from = 10, to = 600, by = 10),
    n_replicates = 100,
    positive_class = "SjD",
    seed = 124
) {
  
  stopifnot(
    is.data.table(data_module),
    length(diag) == nrow(data_module)
  )
  
  set.seed(seed)
  
  diag <- as.factor(diag)
  
  results <- list()
  res_id <- 1
  
  for (batch_size in batch_sizes) {
    message("batch_size : ",batch_size) 
    for (rep in seq_len(n_replicates)) {
      
      idx <- sample(seq_len(nrow(data_module)), size = batch_size, replace = FALSE)
      
      batch_data <- data_module[idx]
      batch_diag <- diag[idx]
      
     
      n_sjd <- sum(batch_diag == "SjD")
      n_ctrl <- sum(batch_diag == "Control")
      
      
      if (n_sjd == 0 || n_ctrl == 0) {
        auc_value <- NA_real_
      } else {
        
        
        ts_dt <- compute_transition_score_from_data(
          data_module = batch_data,
          diag = batch_diag
        )
        
       
        roc_obj <- roc(
          response = batch_diag,
          predictor = ts_dt$TS,
          levels = c("Control", positive_class),
          direction = "<"
        )
        
        auc_value <- as.numeric(auc(roc_obj))
      }
      
      
      results[[res_id]] <- data.table(
        batch_size = batch_size,
        replicate = rep,
        n_SjD = n_sjd,
        n_Control = n_ctrl,
        AUC = auc_value
      )
      
      res_id <- res_id + 1
    }
  }
  
  
  return(rbindlist(results))
}


res <- test_ts_robustness(data_module, diag, batch_sizes = granularity)
saveRDS(res, paste0(folder_for_res,"low_point_nb_test.rds"))
res=readRDS(paste0(folder_for_res,"low_point_nb_test.rds"))

# 4. variance  -----------------------------------------------------------------

res[, AUC_var := var(AUC, na.rm = TRUE), by = batch_size]

#plateau
var_dt <- res[, .(
  AUC_var = var(AUC, na.rm = TRUE)
), by = batch_size][order(batch_size)]

window <- 5  

var_dt[, gain := AUC_var - shift(AUC_var, n = window, type = "lead")]

plateau_start <- var_dt[gain < max(var_dt$AUC_var)/100][1, batch_size]
print(paste("plateau :",  plateau_start))




# 5. plots ----------------------------------------------------------------


plot_var = ggplot(var_dt, aes(batch_size, AUC_var)) +
  geom_point() +
  geom_vline(xintercept = plateau_start, linetype = "dashed", color = "red") +
  theme_bw() +
  labs(
    x = "Batch size",
    y = "AUC variance"
  )+
  scale_x_continuous(breaks=granularity)+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


plot_auc=ggplot(res, aes(x = batch_size, y = AUC)) +
  geom_boxplot(aes(group = batch_size)) +
  geom_jitter(width = 5, alpha = 0.05) +
  geom_vline(xintercept = plateau_start,
             linetype = "dashed", color = "red") +
  theme_bw() +
  labs(
    x = "Batch size",
    y = "AUC"
  ) +
  ylim(min(res$AUC), 1)+
  scale_x_continuous(breaks=granularity)+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


sensivity_analysis_plots = readRDS(paste0(folder_for_res,"sensitivity_analysis_plots.rds"))
a=sensivity_analysis_plots$a
b=sensivity_analysis_plots$b
c=sensivity_analysis_plots$c
d=sensivity_analysis_plots$d


e <- wrap_elements(full = as_grob(plot_auc))
f <- wrap_elements(full = as_grob(plot_var))

comb <- a/b /c/d/e/f +
  plot_layout(ncol = 2, guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(
    plot.tag          = element_text(size = 20, face = "bold"),
    plot.tag.position = c(0, 1)
  )

ggsave(plot = comb, filename=paste0(folder_for_res,"sensitivity_analysis.pdf"),height = 20, width=20)

