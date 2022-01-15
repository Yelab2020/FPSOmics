
#load('data/module3_SKCM_Cli_infor.rda')
#load('data/module3_data_ES_score.rda')
#load('data/module3_mRNAseq.pri.rda')
#source('R/cal.R')
#source('R/FPS_multi_omics_features.R')
#source('R/GIPW_function_omega.R')
#sum.mRNAAll=FPS_multi_omics_features(SKCM_Cli,data_ES,mRNAseq.pri)

#' Calculation of ferroptosis-associated molecular signatures
#'
#' A propensity score matching (PSM) algorithm, which appropriately control the effects of clinical confounding factors, was used to calculate ferroptosis-associated molecular signatures between high and low ferroptosis score group.
#'
#'
#' @param input_Clinical data.frame of the clinical information, the column of barcode is samples and the rest column is each clinical feature(For continuous variables, the name could be age_continuous; For discrete variables, the name could be gender_discrete)
#' @param input_FPS_score ferroptosis score result from FPS function
#' @param input_omics molecular matrix (rownames of the variable must be molecular symbol, each column is a sample)
#' @return data.frame of a propensity score matching (PSM) algorithm result with fdr < 0.05,column of feature.sig is each gene, column of pvalue.sig is p value, column of fdr.sig is fdr, column of coef.sig is coefficient, column of mean0.sig is mean of low group, column of mean1.sig is mean of high group, column of mean0.sig.w is weighted mean of low group, column of mean1.sig.w is weighted mean of high group
#' @export
#'
#' @examples
#' data(m3_input_Clinical,package='FPSOmics')
#' data(m3_input_FPS_score,package='FPSOmics')
#' data(m3_input_mRNA,package='FPSOmics')
#' PSM_result=FPSmultiOmicsSig(m3_input_Clinical,m3_input_FPS_score,m3_input_mRNA)

FPSmultiOmicsSig <- function(input_Clinical,input_FPS_score,input_omics) {

  SKCM_Cli=input_Clinical #clinic infor
  data_ES=input_FPS_score #FPS score
  mRNAseq.pri=input_omics #omisc data
  #scripts.dir=temp_path

  library(dummies)
  library(formula.tools)
  #folder <- "result_PSM"
  #if (!file.exists(folder)) { dir.create(folder) }
  #scripts.dir <- paste(getwd(),"/",folder,sep="")
  #setwd(scripts.dir)

  #这个统一设定
  cancerNames <- c("user") #
  data.type='omics' #

  #==用户输入 1 临床信息
  #SKCM_Cli <- read.table("/work/dy/Collaboration/he/Ferroptosis_website_code/data/module3_clinical_infor_for_deg_input1.tab",header=T,sep='\t')

  #==标记 连续 和 离散 变量---------------需要对混杂因素进行矫正，再计算差异基因
  ContinousFeatures=colnames(SKCM_Cli)[grep('_continuous',colnames(SKCM_Cli))]
  DiscreteFeatures=colnames(SKCM_Cli)[grep('_discrete',colnames(SKCM_Cli))]

  #==用户输入 2 分组信息
  #data_ES <- read.table("/work/dy/Collaboration/he/Ferroptosis_website_code/data/module3_score_for_deg_input2.tab",header=T,sep='\t')

  data_ES=data_ES[,c('barcode','myclusters')]
  SKCM_Cli <- merge(SKCM_Cli,data_ES,by="barcode")
  SKCM_Cli <- SKCM_Cli[!duplicated(SKCM_Cli$barcode),]

  analysis="myclusters" #Oxygen_Content
  sum.mRNAAll <- data.frame()
  ####must exist stum  clinical, stratification, mRNA files.
  cancerNames <- c("user")
  for(cancer in cancerNames){
    data <- SKCM_Cli
    data <- unique(data)
    rownames(data) <- data[,1]
    data <- data[,-1]

    analysis <- "myclusters"
    if(analysis=="myclusters"){
      # convert hypoxic and normoxic to numeric 1,0 to suppress the warning message in lm
      data$myclusters <- ifelse(data$myclusters=="High",1,0)
      colnames(data)[which(colnames(data)=="myclusters")] <- "Z"
    }

    # convert to dummy
    dummy.feature <- setdiff(colnames(data),c("Z",ContinousFeatures))#,"pathologic_stage"))
    if(length(dummy.feature)>0){
      data.dum <- dummy.data.frame(data, names=dummy.feature)
      dummy.list <- attr(data.dum,"dummies")
      rm.col <- c()
      for (i in 1:length(dummy.list))
      {
        rm.col <- c(rm.col, dummy.list[[i]][length(dummy.list[[i]])])
      }
      data.dum <- data.dum[,-rm.col]
      data.dum$X0 <- rep(1, nrow(data.dum))
      #form <- as.formula("Z~.") # should exclude X0
      exclude.col <- match(c("Z","X0"), colnames(data.dum))
      colnames(data.dum) <- gsub(" ", ".", colnames(data.dum))
      form <- as.formula(paste("Z~",paste(colnames(data.dum)[-exclude.col],collapse="+"),sep=""))
    }else{
      data.dum <- data
      data.dum$X0 <- rep(1, nrow(data.dum))
      #form <- as.formula("Z~.") # should exclude X0
      exclude.col <- match(c("Z","X0"), colnames(data.dum))
      colnames(data.dum) <- gsub(" ", ".", colnames(data.dum))
      form <- as.formula(paste("Z~",paste(colnames(data.dum)[-exclude.col],collapse="+"),sep=""))
    }


    # perform calculation
    #source("cal.R")
    library(doMC)
    library(foreach)
    registerDoMC(21)
    # mRNA.exp
    ##==#==用户输入 3  组学谱信息
    #mRNAseq.pri <- read.delim("/work/dy/Collaboration/he/Ferroptosis_website_code/data/module3_mRNAseq.pri_input3.tab",header=T,sep='\t',row.names=1)

    mRNAseq.pri <- t(mRNAseq.pri)
    mRNAseq.pri <- rm.zero.col(mRNAseq.pri)
    folder <- paste(cancer,"_",analysis,sep="")
    if (!file.exists(folder)) { dir.create(folder) }

    mRNAseq.result <- weight.test(data.dum, form, mRNAseq.pri, is.continuous=TRUE,weight=ifelse(analysis=="myclusters","MW","ATT"),mirror.plot=FALSE, cancer, data.type= data.type, outdir='.',perm=FALSE)
    sum.mRNA <- summarize.fdr(mRNAseq.pri, mRNAseq.result, print=TRUE)
    summarize.fdr(mRNAseq.pri, mRNAseq.result, print=TRUE, cutoff=0.05)
    #write.summary(sum.mRNA, cancer, analysis,"mRNA") #可不输出 1
    #write.result(mRNAseq.result, cancer, analysis,"mRNA") #可不输出 2
    #save(mRNAseq.result, file=paste(cancer,"_", analysis,"_result.RData",sep=""))#可不输出 3
    if(length(which(mRNAseq.result$fdr < 0.05)) > 0){
      sum.mRNA <- data.frame(sum.mRNA)
      #sum.mRNA$class <- rep(cancer,times=nrow(sum.mRNA))
      if(nrow(sum.mRNAAll) == 0){
        sum.mRNAAll <- sum.mRNA
      }else{
        sum.mRNAAll <- rbind(sum.mRNAAll,sum.mRNA)
      }
    }
  }

   return(sum.mRNAAll)
#输出表 1
#write.table(sum.mRNAAll,file="mRNAseqlog2.genes.SKCM_output.txt",quote = F,sep="\t",row.names = F)
}
