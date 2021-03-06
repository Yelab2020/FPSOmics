% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FPSmultiOmicsSigEnrich.R
\name{FPSmultiOmicsSigEnrich}
\alias{FPSmultiOmicsSigEnrich}
\title{Function enrichment}
\usage{
FPSmultiOmicsSigEnrich(input_PSM, input_gs_kegg, temp_folder)
}
\arguments{
\item{input_PSM}{data.frame of the PSM result from FPSmultiOmicsSig function}

\item{input_gs_kegg}{the genelist of two column for GSEA enrichemnt, the first column is term and the second column is gene}

\item{temp_folder}{the temporal folder for output result}
}
\value{
table and figures of GO, KEGG and GSEA enrichemnt result
}
\description{
Study the differences of different ferroptosis status in cancer hallmarks, KEGG pathway by Gene Set Enrichment Analysis (GSEA) and GO, KEGG enrichment.
}
\examples{
data(m4_input_gs_kegg,package='FPSOmics')
data(m4_input_PSM,package='FPSOmics')
enrichment_result=FPSmultiOmicsSigEnrich(m4_input_PSM,m4_input_gs_kegg,'temp_folder')

}
