\name{piano-package}
\alias{piano-package}
\alias{piano}
\docType{package}
\title{
  Piano - Platform for Integrative ANalysis of Omics data
}

\description{
 Run gene set analysis with various statistical methods, from different gene level statistics and a wide range of gene-set collections. Furthermore, the Piano package contains functions for combining the results of multiple runs of gene set analyses.
}
\details{
  The Piano package consists of two parts. The major part revolves around gene set analysis (GSA), and the central function for this is \code{\link{runGSA}}. There are some downstream functions (e.g. \code{\link{GSAsummaryTable}} and \code{\link{geneSetSummary}}) that handle the results from the GSA. By running \code{runGSA} multiple times with different settings it is possible to compute consensus gene set scores. Another set of functions (e.g. \code{\link{consensusScores}} and \code{\link{consensusHeatmap}}) take a list of result objects given by \code{runGSA} for this step. The second part of the Piano package contains a set of functions devoted for an easy-to-use approach on microarray analysis (wrapped around the \pkg{affy} and \pkg{\link[limma:limma-package]{limma}} packages), which are constructed to integrate nicely with the downstream GSA part. The starting function in this case is \code{\link{loadMAdata}}.
}

\author{
  Leif Varemo \email{piano.rpkg@gmail.com} and
  Intawat Nookaew \email{piano.rpkg@gmail.com}
}

\seealso{
  \code{\link{runGSA}} and \code{\link{loadMAdata}}
}