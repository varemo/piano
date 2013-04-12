\name{loadGSC}
\alias{loadGSC}
\title{
   Load a gene set collection
}
\description{
   Load a gene set collection, to be used in \code{\link{runGSA}}, in GMT, SBML or SIF format, or optionally from a \code{data.frame}.
}
\usage{
   loadGSC(file, type="auto", addInfo)
}
\arguments{
   \item{file}{
      a character string, giving the name of the file containing the gene set collection. Optionally an object that can be coerced into a two-column data.frame, the first column containing genes and the second gene sets, representing all "gene"-to-"gene set" connections.
   }
   \item{type}{
      a character string giving the file type. Can be either of \code{"gmt"}, \code{"sbml"}, \code{"sif"}. If set to \code{"auto"} the type will be taken from the file extension. If the gene-set collection is loaded into R from another source and stored in a data.frame, it can be loaded with the setting \code{"data.frame"}.
   }
   \item{addInfo}{
      an optional data.frame with two columns, the first containging the gene set names and the second containing additional information for each gene set. Some additional info may load automatically from the different file types.
   }
}
\details{
   This function is used to create a gene-set collection object to be used with \code{\link{runGSA}}. 
   
   The "gmt" files available from the Molecular Signatures Database (\url{http://www.broadinstitute.org/gsea/msigdb/}) can be loaded using \code{loadGSC}. This website is a valuable resource and contains several different collections of gene sets.
   
   By using the functionality of e.g. the \code{biomaRt} package, a gene-set collection with custom gene names (matching the statistics used in \code{\link{runGSA}}) can easily be compiled into a two-column data.frame (column order: genes, gene sets) and loaded with \code{type="data.frame"}. 
   
If a sif-file is used it is assumed that the first column contains gene sets and the third column contains genes. 
   
In order to load an SBML-file it is required that libSBML and \code{rsbml} is installed. Note that the SBML loading is an experimental feature and is highly dependent on the version and format of the SBML file and requires it to contain gene association. By examining the returned \code{GSC} object it is easy to see if the correct gene sets were loaded.
}
\value{
   A list like object of class \code{GSC} containing two elements. The first is \code{gsc}, a list of the gene sets, each element a character vector of genes. The second element is \code{addInfo}, a data.frame containing the optional additional information.
   
}
\author{
   Leif Varemo \email{piano.rpkg@gmail.com} and
  Intawat Nookaew \email{piano.rpkg@gmail.com}
}
\seealso{
   \pkg{\link{piano}}, \code{\link{runGSA}}
}
\examples{
   # Randomly generated gene sets:
   g <- sort(paste("g",floor(runif(100)*500+1),sep=""))
   g <- c(g,sort(paste("g",floor(runif(900)*1000+1),sep="")))
   g <- c(g,sort(paste("g",floor(runif(1000)*2000+1),sep="")))
   s <- paste("s",floor(rbeta(2000,0.9,1.7)*50+1),sep="")
   
   # Make data.frame:
   gsc <- cbind(g,s)
   
   # Load gene set collection from data.frame:
   gsc <- loadGSC(gsc)
}
