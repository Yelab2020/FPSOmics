% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FPS.R
\name{FPS}
\alias{FPS}
\title{Calculate of ferroptosis score}
\usage{
FPS(input_exp_data)
}
\arguments{
\item{input_exp_data}{expression matrix (rownames of the variable must be gene symbol, each column is a sample)}
}
\value{
data.frame of ferroptosis score, column of pro is pro-ferroptosis score, column of anti is anti-ferroptosis score, column of FPS is Ferroptosis score, column of barcode is samples, column of myclusters is high and low group divided by median ferroptosis score
}
\description{
We firstly calculate enrichment score (ES) of ferroptosis-promote genes(pro-FRGs) and ferroptosis-inhibit genes(anti-FRGs) using single sample gene set enrichment analysis (ssGSEA), then ferroptosis score was defined by the ES of pro-FRGs minus anti-FRGs.
}
\examples{
data(m1_input_mRNA,package='FPSOmics')
FPS_score=FPS(m1_input_mRNA)
}
