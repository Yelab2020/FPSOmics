#' Calculate of ferroptosis score
#'
#' We firstly calculate enrichment score (ES) of ferroptosis-promote genes(pro-FRGs) and ferroptosis-inhibit genes(anti-FRGs) using single sample gene set enrichment analysis (ssGSEA), then ferroptosis score was defined by the ES of pro-FRGs minus anti-FRGs.
#'
#' @param input_exp_data  expression matrix (rownames of the variable must be gene symbol, each column is a sample)
#'
#' @return data.frame of ferroptosis score, column of pro is pro-ferroptosis score, column of anti is anti-ferroptosis score, column of FPS is Ferroptosis score, column of barcode is samples, column of myclusters is high and low group divided by median ferroptosis score
#' @export
#'
#' @examples
#' data(m1_input_mRNA,package='FPSOmics')
#' FPS_score=FPS(m1_input_mRNA)

FPS <- function(input_exp_data) {
  library('magrittr')

  ##genelist     input 1
  data(m1_input_genelist,package='FPSOmics')

  ## input 2
  melanoma_m <- input_exp_data

  data_ES <- GSVA::gsva(as.matrix(melanoma_m),genelist,method="ssgsea")
  data_ES <- t(data_ES) %>% data.frame()
  data_ES$FPS <- data_ES$pro-data_ES$anti
  data_ES$barcode=rownames(data_ES)
  data_ES$myclusters <- ifelse(data_ES$FPS>=median(data_ES$FPS),"High","Low")
  #data_ES_zscore <- scale(data_ES) %>% data.frame()

  return(data_ES)

}
