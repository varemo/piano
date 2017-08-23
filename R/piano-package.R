#' Piano - Platform for Integrative ANalysis of Omics data
#' 
#' Run gene set analysis with various statistical methods, from different gene
#' level statistics and a wide range of gene-set collections. Furthermore, the
#' Piano package contains functions for combining the results of multiple runs
#' of gene set analyses.
#' 
#' The Piano package consists of two parts. The major part revolves around gene
#' set analysis (GSA), and the central function for this is
#' \code{\link{runGSA}}. There are some downstream functions (e.g.
#' \code{\link{GSAsummaryTable}} and \code{\link{geneSetSummary}}) that handle
#' the results from the GSA. By running \code{runGSA} multiple times with
#' different settings it is possible to compute consensus gene set scores.
#' Another set of functions (e.g. \code{\link{consensusScores}} and
#' \code{\link{consensusHeatmap}}) take a list of result objects given by
#' \code{runGSA} for this step. The second part of the Piano package contains a
#' set of functions devoted for an easy-to-use approach on microarray analysis
#' (wrapped around the \pkg{affy} and \pkg{\link[limma:limma-package]{limma}}
#' packages), which are constructed to integrate nicely with the downstream GSA
#' part. The starting function in this case is \code{\link{loadMAdata}}.
#' 
#' @name piano-package
#' @aliases piano-package piano
#' @docType package
#' @author Leif Varemo \email{piano.rpkg@@gmail.com} and Intawat Nookaew
#' \email{piano.rpkg@@gmail.com}
#' @seealso \code{\link{runGSA}} and \code{\link{loadMAdata}}
NULL

#' Random input data for gene set analysis
#' 
#' This data set is completely randomly generated and contains p-values for
#' 2000 genes, fold-changes for those genes and a gene set collection giving
#' the connection between genes and 50 gene sets. Only intended to be used as
#' example data for \code{\link{runGSA}}.
#' 
#' 
#' @name gsa_input
#' @docType data
#' @format A list containing 3 elements: gsa_input$pvals and
#' gsa_input$directions are numeric vectors, gsa_input$gsc is a two-column
#' matrix with gene names in the first column and gene set names in the second.
#' @keywords datasets
NULL

#' Gene set analysis result data
#' 
#' This data set contains gene set analysis results, as returned by the
#' \code{\link{runGSA}} function, that is used as example data for downstream
#' functions. The input data to \code{runGSA} was randomly generated and is
#' accessible through \code{data(gsa_input)}.
#' 
#' 
#' @name gsa_results
#' @docType data
#' @format A list where each element is an object of class GSAres, as returned
#' by \code{runGSA}.
#' @keywords datasets
NULL



