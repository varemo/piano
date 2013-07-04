\name{diffExp}
\alias{diffExp}
\title{
  Perform differential expression analysis
}
\description{
  Identifies differentially expressed genes by using the linear model approach of
  \pkg{limma}. Optionally produces a Venn diagram, heatmap and Polar plot.
}
\usage{
diffExp(arrayData, contrasts, chromosomeMapping, fitMethod = "ls",
        adjustMethod = "fdr", significance = 0.001,
        plot = TRUE, heatmapCutoff = 1e-10, colors=c("red","green","blue","yellow",
        "orange","purple","tan","cyan","gray60","black"),
        save = FALSE, verbose = TRUE)
}
\arguments{
  \item{arrayData}{
    an object of class \code{ArrayData}.
}
  \item{contrasts}{
    a character vector giving the contrasts to be tested for differential expression.
    Use \code{\link{extractFactors}} to get allowed contrasts.
}
  \item{chromosomeMapping}{
    character string giving the name of the chromosome mapping file, or an object of
    class \code{data.frame} or similar containing the chromosome mapping. Required for
    the Polar plot if the \code{ArrayData} object lacks annotation information. See
    details below.
}
  \item{fitMethod}{
    character string giving the fitting method used by \code{lmFit}. Can be either
    \code{"ls"} for least squares (default) or \code{"robust"} for robust regression.
}
  \item{adjustMethod}{
    character string giving the method to use for adjustment of multiple testing.
    Can be \code{"holm"}, \code{"hochberg"}, \code{"hommel"}, \code{"bonferroni"},
    \code{"BH"}, \code{"BY"}, \code{"fdr"} (default) or \code{"none"}. See \code{p.adjust}
    for details.
}
  \item{significance}{
    number giving the significance cutoff level for the Venn diagram.
    Defaults to \code{0.001}.
}
  \item{plot}{
    should plots be produced? Set either to \code{TRUE} (default) or \code{FALSE} to control all plots, or to a character vector with any combination of \code{"venn"}, \code{"heatmap"} and \code{"polarplot"}, to control the single plots (e.g. \code{plot=c("venn","polarplot")} or \code{plot="heatmap"}).
}
  \item{heatmapCutoff}{
    number giving the significance cutoff level for the heatmap. Defaults to
    \code{1e-10}.
}
  \item{colors}{
    character vector of colors to be used by the Venn diagram and Polar plot.
}
  \item{save}{
    should the figures and p-values be saved? Defaults to \code{FALSE}.
}
  \item{verbose}{
    verbose? Defaults to \code{TRUE}.
}
}
\details{
  This function uses \pkg{\link[limma]{limma}} to calculate p-values measuring differential
  expression in the given \code{contrasts}. The \code{uniqueFactors} given by
  \code{\link{extractFactors}} can be used to define a contrast vector, where each
  element should be a character string on the form
  \code{"uniqueFactorA - uniqueFactorB"}, note the space surrounding the \code{-}.
  (See the example below and for \code{\link{extractFactors}}.)
  
  If appropriate annotation is missing for the \code{ArrayData} object the user
  can suppply this as \code{chromosomeMapping}. This should be either a \code{data.frame}
  or a tab delimited text file and include the columns \emph{chromosome} with the
  chromosome name and \emph{chromosome location} containing the starting position
  of each gene. A \code{-} sign can be used to denote the antisense strand but this will
  be disregarded while plotting. The rownames should be \emph{probe IDs} or, if
  using a text file, the first column with a column header should contain the
  \emph{probe IDs}.
  
  Note that the \code{fitMethod="robust"} may need longer time to run.
  
  A Venn diagram can be drawn for up to five contrasts (\code{diffExp()} will use \code{vennDiagram}).
  
  The heatmap shows normalized expression values of the genes that pass the
  \code{heatmapCutoff} in at least one contrast.
  
  The Polar plot sorts the genes according to chromosomal location, for each
  chromosome starting with unknown positions followed by increasing number in
  the \emph{chromosome location} column. Genes which do not map to any chromosome
  are listed as U for unknown. The radial lines in the Polar plot are -log10 scaled
  p-values, so that a longer line means a smaller p-value. This gives an overview
  of the magnitude of differential expression for each contrast.

  Typical usages are:
\preformatted{
# Identify significantly changed genes in 'm1' and 'm2' compared to 'wt':
diffExp(arrayData, contrasts=c("m1 - wt", "m2 - wt"))
}
}
\value{
  A \code{list} with elements:

  \item{pValues}{\code{data.frame} containing p-values for each contrast}
  \item{foldChanges}{\code{data.frame} containing log2 fold changes for each contrast}
  \item{resTable}{a \code{list} with an element for each contrast, each being a \code{data.frame} with full result information}
}
\author{
  Leif Varemo \email{piano.rpkg@gmail.com} and
  Intawat Nookaew \email{piano.rpkg@gmail.com}
}
\seealso{
  \pkg{\link{piano}}, \code{\link{loadMAdata}}, \code{\link{extractFactors}},
  \code{\link{polarPlot}}, \code{\link{runGSA}}, \pkg{\link[limma]{limma}},
  \code{\link{venn}}, \code{\link{heatmap.2}}
}
\examples{
  # Get path to example data and setup files:
  dataPath <- system.file("extdata", package="piano")

  # Load normalized data:
  myArrayData <- loadMAdata(datadir=dataPath, dataNorm="norm_data.txt.gz", platform="yeast2")

  # Perform differential expression analysis:
  pfc <- diffExp(myArrayData, contrasts=c("aerobic_Clim - anaerobic_Clim",
                                          "aerobic_Nlim - anaerobic_Nlim"))
                                          
  # Order the genes according to p-values, for aerobic_Clim vs anaerobic_Clim:
  o <- order(pfc$resTable$'aerobic_Clim - anaerobic_Clim'$P.Value)
  
  # Display statistics for the top 10 significant genes:
  pfc$resTable$'aerobic_Clim - anaerobic_Clim'[o[1:10],]

}
\references{
Smyth, G. K. (2005). Limma: linear models for microarray data. In:
  'Bioinformatics and Computational Biology Solutions using R and
  Bioconductor'. R. Gentleman, V. Carey, S. Dudoit, R. Irizarry, W.
  Huber (eds), Springer, New York, pages 397--420.
}


