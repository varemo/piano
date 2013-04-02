\name{extractFactors}
\alias{extractFactors}
\title{
  Extracts \code{ArrayData} factors
}
\description{
  Extracts the factors, given by an \code{ArrayData} object, that can be used by
  \code{\link{diffExp}}
}
\usage{
extractFactors(arrayData)
}
\arguments{
  \item{arrayData}{
    an \code{ArrayData} object.
  }
}
\value{
  A \code{list} with elements:
  
  \item{factors}{Assigns one factor to each array}
  \item{uniqueFactors}{The unique factors that can be used to form contrasts}
}
\author{
  Leif Varemo \email{piano.rpkg@gmail.com} and
  Intawat Nookaew \email{piano.rpkg@gmail.com}
}
\seealso{
  \pkg{\link{piano}}, \code{\link{diffExp}}

}
\examples{
  # Get path to example data and setup files:
  dataPath <- system.file("extdata", package="piano")

  # Load normalized data:
  myArrayData <- loadMAdata(datadir=dataPath, dataNorm="norm_data.txt.gz", platform="yeast2")

  #Extract the factors that can be used in the call to diffExp:
  extractFactors(myArrayData)
}


