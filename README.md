
## 1. Introduction
  FPSOmics is an R package to calculate ferroptosis score (FPS) and multidimensional ferroptosis-associated molecular signatures. 
 FPSOmics designed a novel score model based on ferroptosis-related genes(FRGs)[1] using single sample gene set enrichment analysis (ssGSEA)[2]. This package also include a propensity score matching (PSM) algorithm[3], which appropriately control the effects of clinical confounding factors, to calculate ferroptosis-associated molecular signatures.
  The full FPSOmics documentation is available in the package. To reach the user’s Guide, install the FPSOmics package and load it into an R session by library (FPSOmics). And they can get help by help (FPS), help (FPSsurvival), help (FPSmultiOmicsSig) and help (FPSmultiOmicsSigEnrich) to see documentation of each function. Also description of each function is including in docs/FPSOmics_0.1.0.pdf.

## 2. Installation and requirement
Before installing FPSOmics, please properly install all dependencies
```{r echo = TRUE}
#library(magrittr)
#library(survival)
#library(survminer)
#library(ggplot2)
#library(grid)
#library(gridExtra)
#library(gtable)
#library(org.Hs.eg.db)
#library(KEGG.db)
#library(fgsea)
#library(enrichplot)
#library(clusterProfiler)
#library(magrittr)
#library(dummies)
#library(formula.tools)
#library(doMC)
#library(foreach)
#library(gdata)
#library(gtools)
#library(gmodels)
#library(MASS)
#library(gplots)
#library(formula.tools)
#library(Hmisc) 
#library(Matrix)
#library(dplyr)
```
Installation of FPSOmics

You can install FPSOmics.zip from local path and directly from github.
```{r echo = TRUE}
devtools::install_github('Yelab2020/FPSOmics')
install.packages('FPSOmics.zip", repos = NULL, type = "source")
library(FPSOmics)
```

## 3. Functional modules
### 3.1 Calculate ferroptosis score
We  firstly filtered out to 32 core ferroptosis-related genes (FRGs) from genelist in FerrDb[1].Then, we calculate enrichment score (ES) of  ferroptosis-promote genes (pro-FRGs) and ferroptosis-inhibit genes(anti-FRGs) using single sample gene set enrichment analysis (ssGSEA)[2], then ferroptosis score (FPS) was defined by the ES of pro-FRGs minus anti-FRGs.
```{r}
print('FPS')
data(m1_input_mRNA,package='FPSOmics') #load mRNA matrix
FPS_score=FPSOmics::FPS(m1_input_mRNA)
head(FPS_score)
```

### 3.2 Evaluate prognostic significance
We firstly calculate enrichment score (ES) of ferroptosis-promote genes (pro-FRGs) and ferroptosis-inhibit genes (anti-FRGs) using single sample gene set enrichment analysis (ssGSEA), then ferroptosis score (FPS) was defined by the ES of pro-FRGs minus anti-FRGs. Then Kaplan-Meier survival analyses between high-FPS and low-FPS were also assessed using log rank test.
```{r}
print('FPS-survival')
data(m1_input_mRNA,package='FPSOmics')#load mRNA matrix
data(m2_input_survival_data,package='FPSOmics') #load survival time
#FPS_survival_result=FPSOmics::FPSsurvival(m1_input_mRNA,m2_input_survival_data)

```

### 3.3 Multidimensional ferroptosis-associated molecular signatures
A propensity score matching (PSM) algorithm[3], which appropriately control the effects of clinical confounding factors, was used to calculate ferroptosis-associated molecular signatures between high and low ferroptosis score group.
```{r}
print('FPS-multiOmics')
data(m3_input_Clinical,package='FPSOmics') #load clinical information
data(m3_input_FPS_score,package='FPSOmics') #load ferroptosis score (FPS) result
data(m3_input_mRNA,package='FPSOmics') #load omics data for mRNA
PSM_result=FPSOmics::FPSmultiOmicsSig(m3_input_Clinical,m3_input_FPS_score,m3_input_mRNA)
```


### 3.4 Function enrichment
Study function of ferroptosis-associated molecular signatures and the differences of different ferroptosis status in cancer hallmarks, KEGG pathway by Gene Set Enrichment Analysis (GSEA)[4] and GO, KEGG enrichment by clusterprofile[5].
```{r}
print('FPS-multiOmics-function')
data(m4_input_gs_kegg,package='FPSOmics')#load geneset for GSEA analysis
data(m4_input_PSM,package='FPSOmics')#load result of PSM analysis
#enrichment_result=FPSOmics::FPSmultiOmicsSigEnrich(m4_input_PSM,m4_input_gs_kegg,'temp_folder')

```

## 4. Summary
FPS provides four function modules which allows users to analysis of ferroptosis score and their-related multi-omics molecular signatures. There may still be gaps in the toolkit, but we'll work on more features later. If you have any suggestions or questions, please feel free to contact us (yihe0902@163.com or yudong123@sjtu.edu.cn).

## 5. Session information
```{r}
sessionInfo()
```
## 6. Citing FPSOmics
If FPSOmics R package is used in your research, please cite: Multi-omics characterization and therapeutic liability of ferroptosis in melanoma

## 7. REFERENCES
[1]. Zhou, N. & Bao, J. FerrDb: a manually curated resource for regulators and markers of ferroptosis and ferroptosis-disease associations. Database 2020, 21 (2020).

[2]. Hänzelmann, S., Castelo, R. & Guinney, J. GSVA: gene set variation analysis for microarray and RNA-Seq data. BMC Bioinformatics 14, 7 (2013)

[3]. Ye, Y. et al. Characterization of hypoxia-associated molecular features to aid hypoxia-targeted therapy. Nat. Metab. 1, 431–444 (2019).

[4]. Subramanian, A. et al. Gene set enrichment analysis: A knowledge-based approach for interpreting genome-wide expression profiles. Proc. Natl. Acad. Sci. U. S. A. 102, 15545–15550 (2005)

[5]. Yu, G., Wang, L.-G., Han, Y. & He, Q.-Y. clusterProfiler: an R package for comparing biological themes among gene clusters. OMICS 16, 284–287 (2012).
