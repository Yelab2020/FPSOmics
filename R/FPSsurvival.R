
#' Calculation of ferroptosis score and survival analysis between high and low ferroptosis score groups.
#'
#' We firstly calculate enrichment score (ES) of ferroptosis-promote genes(pro-FRGs) and ferroptosis-inhibit genes(anti-FRGs) using single sample gene set enrichment analysis (ssGSEA), then ferroptosis score was defined by the ES of pro-FRGs minus anti-FRGs. Secondly, survival analysis between high and low ferroptosis score groups.
#'
#' @param input_exp_data expression matrix (rownames of the variable must be gene symbol, each column is a sample)
#' @param input_survival_data  data.frame of the survival time information, the name of first colume is barcode which means samples, the name of second colume is OS which means censos status, the name of third colume is OS.time which means survival time
#'
#' @return data.frame of ferroptosis score, column of pro is pro-ferroptosis score, column of anti is anti-ferroptosis score, column of FPS is Ferroptosis score, column of barcode is samples, column of myclusters is high and low group divided by median ferroptosis score; score_MK_plot is MK curve for high and low ferroptosis score
#' @export
#'
#' @examples data(m1_input_mRNA,package='FPSOmics') #load mRNA matrix
#' data(m2_input_survival_data,package='FPSOmics') #load survival time
#' FPS_survival_result=FPSsurvival(m1_input_mRNA,m2_input_survival_data)
#'
FPSsurvival <- function(input_exp_data,input_survival_data) {

  library('magrittr')

  ##genelist     input 1
  data(m1_input_genelist,package='FPSOmics')

  ##    input 2
  input_exp_data=input_exp_data

  data_ES <- GSVA::gsva(as.matrix(input_exp_data),genelist,method="ssgsea")
  data_ES <- t(data_ES) %>% data.frame()
  data_ES$FPS <- data_ES$pro-data_ES$anti
  data_ES$barcode=rownames(data_ES)
  data_ES$myclusters <- ifelse(data_ES$FPS>=median(data_ES$FPS),"High","Low")

  print('score done')
  ##===========================================================================================================================================================================================
  ## input 3
  survival_data <- input_survival_data

  melcli_group <- merge(data_ES,survival_data,by="barcode")
  melcli_group$Group <- ifelse(melcli_group$FPS>=median(melcli_group$FPS),"high FPS","low FPS")

  ## 绘制生存曲线
  library(survival)
  library(survminer)

  diff=survdiff(Surv(OS.time, OS) ~ Group,data = melcli_group)
  pValue=1-pchisq(diff$chisq,df=1)
  pValue=signif(pValue,2)
  pValue=format(pValue, scientific = TRUE)

  melcli_group$Group <- factor(melcli_group$Group,levels = c("low FPS", "high FPS"))
  cox.res=coxph(Surv(OS.time, OS)~ Group,data = melcli_group)
  HR <- signif(exp(coef(cox.res)),2)
  CI <- signif(exp(confint(cox.res)),2)
  CI_low <- signif(CI[1],2)
  CI_high <- signif(CI[2],2)

  fit <- survfit(Surv(OS.time, OS) ~ Group, data = melcli_group)


  #pdf(file="/work/dy/Collaboration/he/Ferroptosis_website_code/result/module2_outputsurvival.pdf",onefile = FALSE,
  #    width = 5,
  #    height =5)
  score_MK_plot=ggsurvplot(fit,
                           data=melcli_group,
                           conf.int=F,
                           pval=paste0("p=",pValue,"\nHR=",paste0(HR,"(",CI_low,"-",CI_high,")")),
                           pval.size=4,
                           risk.table=TRUE,
                           legend.labs=c("FPS-low","FPS-high"), #, "Cluster4"
                           legend.title="",
                           xlab="Time(years)",
                           ylab="Overall Survival",
                           break.time.by = 5,
                           risk.table.title="",
                           palette=c("lightblue","red"), #, "darkblue"
                           risk.table.height=.3)
  #dev.off()
  print('survival done')
  return(list(data_ES=data_ES,score_MK_plot=score_MK_plot))

}
