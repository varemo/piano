\name{loadMAdata}
\alias{loadMAdata}
\title{
  Load and preprocess microarray data
}
\description{
  Loads, preprocesses and annotates microarray data to be further used by downstream
  functions in the \pkg{\link{piano}} package.
}
\usage{
  loadMAdata(datadir = getwd(), setup = "setup.txt", dataNorm,
           platform = "NULL", annotation, normalization = "plier",
           filter = TRUE, verbose = TRUE, \dots)
}
\arguments{
  \item{datadir}{
    character string giving the directory in which to look for the data. Defaults
    to \code{getwd()}.
}
  \item{setup}{
    character string giving the name of the file containing the experimental setup,
    or an object of class \code{data.frame} or similar containing the experimental
    setup. Defaults to \code{"setup.txt"}, see details below for more information.
}
  \item{dataNorm}{
    character string giving the name of the normalized data, or an object of class
    \code{data.frame} or similar containing the normalized data. Only to be used if
    the user wishes to start with normalized data rather then CEL files.
}
  \item{platform}{
    character string giving the name of the platform, can be either \code{"yeast2"}
    or \code{NULL}. See details below for more information.
}
  \item{annotation}{
    character string giving the name of the annotation file, or an object of class
    \code{data.frame} or similar containing the annotation information. The annotation
    should consist of the columns \emph{Gene name}, \emph{Chromosome} and
    \emph{Chromosome location}. Not required if \code{platform="yeast2"}.
}
  \item{normalization}{
    character string giving the normalization method, can be either \code{"plier"},
    \code{"rma"} or \code{"mas5"}. Defaults to \code{"plier"}.
}
  \item{filter}{
    should the data be filtered? If \code{TRUE} then probes not present in the
    annotation will be discarded. Defaults to \code{TRUE}.
}
  \item{verbose}{
    verbose? Defaults to \code{TRUE}.
}
  \item{\dots}{
    additional arguments to be passed to \code{ReadAffy}.
}
}
\details{
  This function requires at least two inputs: (1) data, either CEL files in the
  directory specified by \code{datadir} or normalized data specified by \code{dataNorm},
  and (2) experimental setup specified by \code{setup}.

  The setup shold be either a tab delimited text file with column headers or a
  \code{data.frame}. The first column should contain the names of the CEL files
  or the column names used for the normalized data, please be sure to use names valid as 
  column names, e.g. avoid names starting with numbers. Additional columns should
  assign attributes in some category to each array. (For an example run the example
  below and look at the object \code{myArrayData$setup}.)
  
  The \pkg{piano} package is customized for yeast 2.0 arrays and annotation will
  work automatically, if the cdfName of the arrays equals \emph{Yeast_2}. If using
  normalized yeast 2.0 data as input, the user needs to set the argument
  \code{platform="yeast2"} to tell the function to use yeast annotation. If other
  platforms than yeast 2.0 is used, set \code{platform=NULL} (default) and supply
  appropriate annotation by the argument \code{annotation}. Note that the cdfName
  will override \code{platform}, so it can still be set to \code{NULL} for yeast
  2.0 CEL files. Note also that \code{annotation} overrides \code{platform}, so
  if the user wants to use an alternative annotation for yeast, this can be done
  simply by specifying this in \code{annotation}.
  
  The annotation should have the column headers \emph{Gene name}, \emph{Chromosome} and
  \emph{Chromosome location}. The \emph{Gene name} is used in the heatmap in
  \code{diffExp} and the \emph{Chromosome} and \emph{Chromosome location} is used
  by the \code{polarPlot}. The rownames (or first column if using a text file)
  should contain the \emph{probe IDs}. If using a text file the first column should
  have the header \emph{probeID} or similar. The filtering step discards all probes
  not listed in the annotation.
  
  Normalization is performed on all CEL file data using one of the Affymetrix methods:
  PLIER (\code{"plier"}) as implemented by \code{\link[plier:justPlier]{justPlier}},
  RMA (Robust Multi-Array Average) (\code{"rma"}) expression measure as implemented
  by \code{\link[affy:rma]{rma}} or
  MAS 5.0 expression measure \code{"mas5"} as implemented by \code{\link[affy:mas5]{mas5}}.
  
  It is possible to pass additional arguments to \code{\link[affy:read.affybatch]{ReadAffy}}, e.g.
  \code{cdfname} as this might be required for some types of CEL files.

}
\value{
  An \code{ArrayData} object (which is essentially a \code{list}) with the
  following elements:

  \item{dataRaw}{raw data as an AffyBatch object}
  \item{dataNorm}{\code{data.frame} containing normalized expression values}
  \item{setup}{\code{data.frame} containing experimental setup}
  \item{annotation}{\code{data.frame} containing experimental setup}
  
  Depending on input arguments the \code{ArrayData} object may not include
  \code{dataRaw} and/or \code{annotation}.

}
\author{
  Leif Varemo \email{piano.rpkg@gmail.com} and
  Intawat Nookaew \email{piano.rpkg@gmail.com}
}
\seealso{
  \pkg{\link{piano}}, \code{\link{runQC}}, \code{\link{diffExp}},
  \code{\link[affy:read.affybatch]{ReadAffy}}, \code{\link[affy:expresso]{expresso}}, \code{\link[plier:justPlier]{justPlier}},
  \code{\link[yeast2.db]{yeast2.db}}
}
\examples{
  # Get path to example data and setup files:
  dataPath <- system.file("extdata", package="piano")

  # Load normalized data:
  myArrayData <- loadMAdata(datadir=dataPath, dataNorm="norm_data.txt.gz", platform="yeast2")

  # Print to look at details:
  myArrayData

}
\references{
Gautier, L., Cope, L., Bolstad, B. M., and Irizarry, R. A.
  affy - analysis of Affymetrix GeneChip data at the probe level.
  \emph{Bioinformatics.} \bold{20}, 3, 307-315 (2004).
}


