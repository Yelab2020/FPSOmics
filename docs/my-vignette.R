## ----echo = TRUE--------------------------------------------------------------
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

## ----echo = TRUE--------------------------------------------------------------
#devtools::install_github('Yelab2020/FPSOmics')
#install.packages('FPSOmics.tar.gz", repos = NULL, type = "source")
#library(FPSOmics)

## -----------------------------------------------------------------------------
print('FPS')
#data(m1_input_mRNA,package='FPSOmics')
#FPS_score=FPSOmics::FPS(m1_input_mRNA)
#head(FPS_score)

## -----------------------------------------------------------------------------
print('FPS-survival')
#data(m1_input_mRNA,package='FPSOmics')#load mRNA matrix
#data(m2_input_survival_data,package='FPSOmics') #load survival time
#FPS_survival_result=FPSOmics::FPSsurvival(m1_input_mRNA,m2_input_survival_data)


## -----------------------------------------------------------------------------
print('FPS-multiOmics')
#data(m3_input_Clinical,package='FPSOmics')
#data(m3_input_FPS_score,package='FPSOmics')
#data(m3_input_mRNA,package='FPSOmics')
#PSM_result=FPSOmics::FPSmultiOmicsSig(m3_input_Clinical,m3_input_FPS_score,m3_input_mRNA)

## -----------------------------------------------------------------------------
print('FPS-multiOmics-function')
#data(m4_input_gs_kegg,package='FPSOmics')
#data(m4_input_PSM,package='FPSOmics')
#enrichment_result=FPSOmics::FPSmultiOmicsSigEnrich(m4_input_PSM,m4_input_gs_kegg,'temp_folder')


## -----------------------------------------------------------------------------
sessionInfo()

