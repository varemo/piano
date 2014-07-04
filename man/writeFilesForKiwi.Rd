\name{writeFilesForKiwi}
\alias{writeFilesForKiwi}

\title{
   Write files for Kiwi
}
\description{
   Given a single object or a list of objects of class GSAres, extract the information needed for visualization in the external python function Kiwi and write it to files that can be used as input.
}
\usage{
   writeFilesForKiwi(gsaRes, label="", overwrite=FALSE)
                    }

\arguments{
  \item{gsaRes}{
      either an object of class \code{GSAres} or a list where each element is an object of class \code{GSAres}, as returned by the \code{runGSA} function.
}
\item{label}{
      a character string that will be appended to the names of the resulting files.
}
\item{overwrite}{
      a logical, whether or not to overwrite existing files with identical names.
}
}
\details{
  This function takes the result from a gene set analysis as returned by the \code{\link{runGSA}} function and writes three files that can be directly used as input to Kiwi. Kiwi is a external function i python that can be used for network-based visualization of the GSA results (http://sysbio.se/kiwi).
}
\value{
   Three files are written in the current directory. GSC.txt contains the gene-gene set associations, i.e. the gene set collection. GLS.txt contains the gene-level statistics. GSS.txt contains the gene set statistics.}

\author{
   Leif Varemo \email{piano.rpkg@gmail.com}
}

\seealso{
   \pkg{\link{piano}}, \code{\link{runGSA}}, \code{\link{networkPlot}}
}
\examples{

   # Load some example GSA results:
   data(gsa_results)
   
   # Write the files:
   writeFilesForKiwi(gsa_results,"exp1")
}


