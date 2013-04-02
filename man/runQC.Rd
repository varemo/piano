\name{runQC}
\alias{runQC}
\title{
  Run quality control
}
\description{
  Performs a set of quality control methods and produces the results as figures.
}
\usage{
runQC(arrayData, rnaDeg = TRUE, nuseRle = TRUE, hist = TRUE,
      boxplot = TRUE, pca = TRUE, colorFactor = 1,
      colors = c("red", "green", "blue", "yellow", "orange",
      "purple", "tan", "cyan", "gray60", "black", "white"),
      save = FALSE, verbose = TRUE)
}
\arguments{
  \item{arrayData}{
    an object of class \code{ArrayData}.
}
  \item{rnaDeg}{
    should RNA degradation be detected? Defaults to \code{TRUE}.
}
  \item{nuseRle}{
    should Normalized Unscaled Standard Errors (NUSE) and Relative Log Expressions
    (RLE) be calculated? Defaults to \code{TRUE}.
}
  \item{hist}{
    produce histograms of expression values? Defaults to \code{TRUE}.
}
  \item{boxplot}{
    produce boxplots of expression values? Defaults to \code{TRUE}.
}
  \item{pca}{
    should PCA be run? Defaults to \code{TRUE}.
}
  \item{colorFactor}{
    a number specifying which column of the setup (given by the \code{ArrayData}
    object) should be used for coloring information for the PCA. Defaults to 1.
}
  \item{colors}{
    a character vector of colors to be used in the PCA plot.
}
  \item{save}{
    should the figures be saved? Defaults to \code{FALSE}.
}
  \item{verbose}{
    verbose? Defaults to \code{TRUE}.
}
}
\details{
  This function is essentially a wrapper for various available quality control
  functions for \code{AffyBatch} objects and normalized microarray data. RNA
  degradation (\code{rnaDeg=TRUE}) and NUSE & RLE (\code{nuseRle=TRUE}) require
  raw data (a \code{dataRaw} element in the \code{ArrayData} object).
  
  The PCA uses \code{\link{prcomp}} on centralized normalized data.

  Typical usages are:
\preformatted{
# Run all quality controls:
runQC(arrayData)
}
}
\author{
  Leif Varemo \email{piano.rpkg@gmail.com} and
  Intawat Nookaew \email{piano.rpkg@gmail.com}
}
\seealso{
  \pkg{\link{piano}}, \code{\link{loadMAdata}}, \code{\link{diffExp}},
  \code{\link[affy:AffyRNAdeg]{AffyRNAdeg}}, \code{\link[affyPLM:fitPLM]{fitPLM}}, \code{\link[affy:AffyBatch-class]{AffyBatch}},
  \code{\link{prcomp}}
}
\examples{
  # Get path to example data and setup files:
  dataPath <- system.file("extdata", package="piano")

  # Load normalized data:
  myArrayData <- loadMAdata(datadir=dataPath, dataNorm="norm_data.txt.gz", platform="yeast2")

  # Run PCA only:
  runQC(myArrayData,rnaDeg=FALSE, nuseRle=FALSE, hist=FALSE, boxplot=FALSE)
}
\references{
Brettschneider J, Collin F, Bolstad BM, and Speed TP. Quality assessment for
short oligonucleotide arrays. \emph{Technometrics.} (2007) In press
}


