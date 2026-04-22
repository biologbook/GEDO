folder_script_to_source = "/home/clem/GEDO/R/"
source(file = paste0(folder_script_to_source,"GEDO.R"))
source(file = paste0(folder_script_to_source,"functions_article.R"))

detach_all_packages()
packages_list = c("data.table", "tryCatchLog","rlist",
                  "dbscan", "progressr",
                  "pROC", 
                  "devtools","progress","ggplot2","patchwork","cowplot","rstatix")

install_and_load(packages_list)

folder_for_res = "/home/clem/GEDO/results/"
folder_to_data ="/home/clem/GEDO/results/"


res=readRDS(paste0(folder_for_res, "module_auc_comparison_by_config.rds"))


# 2. Data.tables -------------------------------------------------------------


dt_long <- melt(
  res,
  id.vars = "pathway",
  variable.name = "config",
  value.name = "auc"
)

dt_long[, method := ifelse(grepl("^GEDOcorr", config), "GEDOcorr", "GEDO")]

dt_long[, core_pct := as.numeric(sub(".*_core([0-9.]+)_.*", "\\1", config))]
dt_long[, k_lof    := as.integer(sub(".*_klof([0-9]+)_.*", "\\1", config))]
dt_long[, k_knn    := as.integer(sub(".*_kknn([0-9]+)$", "\\1", config))]

#Vérification : 
dt_long[, .N, by = .(method, core_pct, k_lof, k_knn)]


#Vérif : 
dt_long[
  ,
  .(
    mean_auc = mean(auc, na.rm = TRUE),
    sd_auc   = sd(auc, na.rm = TRUE)
  ),
  by = .(method, core_pct, k_lof, k_knn)
][order(method, core_pct, k_lof, k_knn)]


#------------K KNN

dt_box_kknn <- dt_long[k_lof==30 & core_pct==0.9,.(auc, method, k_knn,pathway)]
dt_box_kknn[, k_knn := factor(k_knn, levels = c(5, 10, 20, 30))]
dt_box_kknn[, method := factor(method, levels = c("GEDO", "GEDOcorr"))]



#------------K lof

dt_box_klof <- dt_long[k_knn==10 & core_pct==0.9,.(auc, method, k_lof,pathway)]
dt_box_klof[, k_lof := factor(k_lof, levels = c(5, 10, 20, 30))]
dt_box_klof[, method := factor(method, levels = c("GEDO", "GEDOcorr"))]


#---------CORE PCT
dt_box_core <- dt_long[k_knn==10 & k_lof==30,.(auc, method, core_pct,pathway)]
dt_box_core[, method := factor(method, levels = c("GEDO", "GEDOcorr"))]
dt_box_core[, core_pct := factor(core_pct,
                                 levels = c(1.0, 0.9, 0.8, 0.7),
                                 labels = c("100", "90", "80", "70"))]





# 4. ANOVA + plot ---------------------------------------------------------

gedo_cols <- c(
  "GEDO"     = "#0D0887FF",
  "GEDOcorr" = "#E16462FF"
)


fmt_p <- function(p) {
  if (is.na(p)) "NA"
  else if (p < 1e-3) "< 1e-3"
  else sprintf("= %.3g", p)
}

#' Run repeated mesure ANOVA and create a title
#' @param dt data to test
#' @param factor_col groups to compare
#' @param factor_label name for plots
anova_title_repeated <- function(dt, factor_col, factor_label) {
  

  extract_rm_stats <- function(data, factor_col) {
    
    res <- anova_test(
      data = data,
      dv = auc,
      wid = pathway,
      within = !!rlang::sym(factor_col)
    )
    
    tab <- get_anova_table(res)
    
    list(
      r2 = tab$ges[1],   # generalized eta squared (RM effect size)
      p  = tab$p[1]
    )
  }
  
  ## ---------- GEDO ----------
  stats_gedo <- extract_rm_stats(
    data = dt[method == "GEDO"],
    factor_col = factor_col
  )
  
  ## ---------- GEDOcorr ----------
  stats_corr <- extract_rm_stats(
    data = dt[method == "GEDOcorr"],
    factor_col = factor_col
  )
  
  ## ---------- Composite title ----------
  sprintf(
    "%s\nGEDO : ges = %.3f, p %s | GEDOcorr : ges = %.3f, p %s",  #GES = proportion of explained variance (generalized eta squared)
    factor_label,
    stats_gedo$r2, fmt_p(stats_gedo$p),
    stats_corr$r2, fmt_p(stats_corr$p)
  )
}




#----------K-knn
cat("repeated ANOVA K KNN \n")
title_kknn <- anova_title_repeated(
  dt_box_kknn,
  factor_col   = "k_knn",
  factor_label = "Effect of k-KNN on AUC"
)



p_kknn=ggplot(dt_box_kknn, aes(x = k_knn, y = auc, fill = method)) +
  geom_jitter(
    aes(color = method),
    position = position_jitterdodge(
      jitter.width = 0.5,
      dodge.width = 0.8
    ),
    size = 0.6,
    alpha = 0.05
  ) +
  geom_boxplot(
    position = position_dodge(width = 0.8),
    outlier.shape = NA,
    alpha = 0.6
  ) +
  
  scale_fill_manual(values=gedo_cols)+
  scale_color_manual(values=gedo_cols)+
  labs(x = "k-KNN", y = "AUC", title=title_kknn, color="Method", fill="Method") +
  theme_bw()+
  ylim(0.7, 0.975)


#-------------k - lof

cat("repeated ANOVA K LOF \n")

title_klof <- anova_title_repeated(
  dt_box_klof,
  factor_col   = "k_lof",
  factor_label = "Effect of k-LOF on AUC"
)

p_klof=ggplot(dt_box_klof, aes(x = k_lof, y = auc, fill = method)) +
  geom_jitter(
    aes(color = method),
    position = position_jitterdodge(
      jitter.width = 0.5,
      dodge.width = 0.8
    ),
    size = 0.6,
    alpha = 0.05
  ) +
  geom_boxplot(
    position = position_dodge(width = 0.8),
    outlier.shape = NA,
    alpha = 0.6
  ) +
  
  scale_fill_manual(values=gedo_cols)+
  scale_color_manual(values=gedo_cols)+
  labs(x = "k-LOF", y = "AUC", title=title_klof, color="Method", fill="Method") +
  theme_bw() 


#--------------core % 

cat("repeated ANOVA CORE % \n")

title_core <- anova_title_repeated(
  dt_box_core,
  factor_col   = "core_pct",
  factor_label = "Effect of CORE % on AUC"
)

p_core=ggplot(dt_box_core, aes(x = core_pct, y = auc, fill = method)) +
  geom_jitter(
    aes(color = method),
    position = position_jitterdodge(
      jitter.width = 0.5,
      dodge.width = 0.8
    ),
    size = 0.6,
    alpha = 0.05
  ) +
  geom_boxplot(
    position = position_dodge(width = 0.8),
    outlier.shape = NA,
    alpha = 0.6
  ) +
  
  scale_fill_manual(values=gedo_cols)+
  scale_color_manual(values=gedo_cols)+
  labs(x = "CORE %", y = "AUC", title=title_core, color="Method", fill="Method") +
  theme_bw() 
  


a <- wrap_elements(full = as_grob(p_kknn))
b <- wrap_elements(full = as_grob(p_klof))
c <- wrap_elements(full = as_grob(p_core))


comb <- a / b / c   +
  plot_layout(ncol = 3, guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(
    plot.tag          = element_text(size = 20, face = "bold"),
    plot.tag.position = c(0, 1)
  )


comb
ggsave(plot = comb, filename = paste0(folder_for_res,"config_auc_comparisons.png"), height = 8, width=25)




# 5. UMAP GEDO ------------------------------------------------------------

folder_for_data_umap = paste0(folder_for_res,"/grid_search_umap_gedo/")
res_umap=compute_auc_modules_grid_search(folder_for_data = folder_for_data_umap)
saveRDS(res_umap, paste0(folder_for_res, "module_auc_comparison_by_config_umap.rds"))
# res_umap=readRDS("/shared/projects/toposads/finalresult/7_gedo_reviewing/grid_evaluation/module_auc_comparison_by_config_umap.rds")

dt_long_umap <- melt(
  res_umap,
  id.vars = "pathway",
  variable.name = "config",
  value.name = "auc"
)



dt_long_umap[, core_pct := as.numeric(sub(".*_core([0-9.]+)_.*", "\\1", config))]
dt_long_umap[, k_lof    := as.integer(sub(".*_klof([0-9]+)_.*", "\\1", config))]
dt_long_umap[, k_knn :=as.factor(sub(".*_kknn([0-9]+)_.*", "\\1", config))]
dt_long_umap[, umap_ncomp :=as.factor(sub(".*_ncomp([0-9]+).*", "\\1", config))]
dt_long_umap[, k_knn:=factor(k_knn, levels=c("5","10", "20","30"))]
dt_long_umap[, umap_ncomp:=factor(umap_ncomp, levels=c("2","10", "20","30"))]


#' Run ANOVA for UMAP-GEDO
#' @param dt data to test
#' @param factor_test groups to compare
#' @param factor_fixed other parameter fixed
#' @param fixed_value value 
extract_rm_stats_conditional <- function(
    dt,
    factor_test,
    factor_fixed,
    fixed_value
) {
  library(rstatix)
  library(rlang)
  
  dt_sub <- dt[get(factor_fixed) == fixed_value]
  
  res <- anova_test(
    data   = dt_sub,
    dv     = auc,
    wid    = pathway,
    within = !!sym(factor_test)
  )
  
  tab <- get_anova_table(res)
  
  list(
    ges = tab$ges[1],
    p   = tab$p[1],
    n_pathways = length(unique(dt_sub$pathway))
  )
}


#' Run repeated mesure ANOVA and create a title
#' @param dt data to test
#' @param kknn_ref number of neighbors for knn
#' @param ncomp_ref number of components for UMAP
#' @param label title for plots
anova_title_umap_repeated <- function(
    dt,
    kknn_ref = 10,
    ncomp_ref = 10,
    label = "UMAP sensitivity analysis"
) {
  
  stats_kknn <- extract_rm_stats_conditional(
    dt = dt,
    factor_test  = "k_knn",
    factor_fixed = "umap_ncomp",
    fixed_value  = ncomp_ref
  )
  
  stats_ncomp <- extract_rm_stats_conditional(
    dt = dt,
    factor_test  = "umap_ncomp",
    factor_fixed = "k_knn",
    fixed_value  = kknn_ref
  )
  
  sprintf(
    "%s\nk-KNN (ncomp=%s): ges = %.3f, p %s | UMAP N COMP (k-KNN=%s) : ges = %.3f, p %s",
    label,
    ncomp_ref,
    stats_kknn$ges,  fmt_p(stats_kknn$p),
    kknn_ref,
    stats_ncomp$ges, fmt_p(stats_ncomp$p)
  )
}




cat("repeated ANOVA UMAP N COMP \n")

title_umap <- anova_title_umap_repeated(dt = dt_long_umap,label = "Effect of UMAP N COMP and k-KNN on AUC")



p_umap=ggplot(dt_long_umap, aes(x = umap_ncomp, y = auc, fill = k_knn)) +
  geom_jitter(
    aes(color = k_knn),
    position = position_jitterdodge(
      jitter.width = 0.5,
      dodge.width = 0.8
    ),
    size = 0.6,
    alpha = 0.05
  ) +
  geom_boxplot(
    position = position_dodge(width = 0.8),
    outlier.shape = NA,
    alpha = 0.6
  ) +
  labs(x = "UMAP N COMP", y = "AUC", title=title_umap, fill="K KNN", color="K KNN") +
  theme_bw() 
p_umap

d <- wrap_elements(full = as_grob(p_umap))

saveRDS(object = list(a=a,b=b,c=c,d=d), file = paste0(folder_for_res,"sensitivity_analysis_plots.rds"))

comb <- a / b / c /d   +
  plot_layout(ncol = 2, guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(
    plot.tag          = element_text(size = 20, face = "bold"),
    plot.tag.position = c(0, 1)
  )


comb
ggsave(plot = comb, filename = paste0(folder_for_res,"config_auc_comparisons_with_umap.png"), height = 15, width=25)






# 6. Export tables --------------------------------------------------------

cat("Formatting tables  \n")


library(data.table)

summary_gedo <- dt_long[
  ,
  .(
    auc_median = median(auc, na.rm = TRUE),
    auc_sd     = sd(auc, na.rm = TRUE),
    n_pathways = .N
  ),
  by = .(method, core_pct, k_lof, k_knn)
]


summary_gedo[, auc_summary :=
               sprintf("%.3f ± %.3f", auc_median, auc_sd)
]

summary_gedo_table <- summary_gedo[
  ,
  .(
    Method      = method,
    Core_pct    = core_pct,
    K_LOF       = k_lof,
    K_KNN       = k_knn,
    `AUC (median ± sd)` = auc_summary
  )
][order(Method, Core_pct, K_LOF, K_KNN)]



summary_umap <- dt_long_umap[
  ,
  .(
    auc_median = median(auc, na.rm = TRUE),
    auc_sd     = sd(auc, na.rm = TRUE),
    n_pathways = .N
  ),
  by = .(core_pct, k_lof, k_knn, umap_ncomp)
]

summary_umap[, auc_summary :=
               sprintf("%.3f ± %.3f", auc_median, auc_sd)
]

summary_umap_table <- summary_umap[
  ,
  .(
    Method      = "UMAP_GEDO",
    Core_pct    = core_pct,
    K_LOF       = k_lof,
    K_KNN       = k_knn,
    UMAP_ncomp  = umap_ncomp,
    `AUC (median ± sd)` = auc_summary
  )
][order(UMAP_ncomp, K_KNN)]


summary_gedo_table[, UMAP_ncomp := NA_integer_]

final_auc_table <- rbind(
  summary_gedo_table,
  summary_umap_table,
  fill = TRUE
)

setcolorder(
  final_auc_table,
  c("Method", "Core_pct", "K_LOF", "K_KNN", "UMAP_ncomp", "AUC (median ± sd)")
)
setnames(final_auc_table, old=c("K_KNN","K_LOF", "UMAP_ncomp","Core_pct"), new=c("k-KNN", "k-LOF","UMAP N COMP","CORE %"))

write.csv2(final_auc_table, file = paste0(folder_for_res,"grid_search_auc_table.csv"))
