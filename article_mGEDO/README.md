Alpha-Interferons and CD40: complementary pathways revealed by mGEDO in distinct transcripto-methylomic profiles of patients with Sjögren’s disease.
doi : 

# 📘 Abstract

Objectives: Sjögren’s disease (SjD) is clinically and molecularly heterogeneous, which limits robust patient stratification and may contribute to variable responses to targeted therapies. The objective of this study is to quantify alpha-Interferon, CD40 and Interleukin 2 pathways in different SjD subgroups.

Methods: We developed mGEDO, a multi-omic extension of GEDO, to integrate whole-blood RNA-seq and DNA methylation data at the gene-set level using BloodGen3 modules. We applied mGEDO to 702 PRECISESADS participants to get two-dimensional embeddings. We then designed weighted alpha-Interferon-, CD40-, and IL2-related signatures using an automated gene selection method, and computed patient-level scores using z-score summaries and customized regularized canonical correlation analysis.

Results: mGEDO captured a continuous SjD–control molecular landscape and revealed substantial within-group heterogeneity beyond discrete clustering. Alpha-Interferon, CD40, and IL-2 signatures showed distinct but partially overlapping patterns. Alpha-Interferon remained the dominant axis associated with severe biological activity, while CD40 and IL-2 provided complementary information, particularly in IFN-negative subsets. Stratification by combined IFN/CD40 status identified biologically meaningful patient classes, with double-positive patients displaying the strongest immune dysregulation.

Conclusions: mGEDO provides an interpretable framework for integrating transcriptomic and methylomic data and for mapping SjD as a molecular continuum. The complementarity between IFN and CD40/IL2 signatures supports partially uncoupled innate and adaptive immune programs and may improve patient stratification for precision therapeutic targeting.


# 📂 Repository structure
```
GEDO/
├── README.md
├── LICENSE
├── GEDO.Rproj
│
├── R/
│   ├── GEDO/
│   └── mGEDO/
│       └── mGEDO_functions.R
│           Functions implementing the mGEDO method (extension of GEDO),
│           including modified algorithms and additional analysis utilities.
│
├── article_GEDO/
│
├── article_mGEDO/
│   ├── README.md
│   │   Description and instructions to reproduce results for the mGEDO article.
│   │
│   └── scripts/
│       └── mGEDO_compute.qmd
│           Quarto workflow to run the mGEDO method and generate associated results.
│            
```

# 📄 Tutorial 

[▶️ See explicative video](https://github.com/biologbook/GEDO/releases/tag/uploading_GEDO_manim_animation/GEDO_manim_animation.mp4)

