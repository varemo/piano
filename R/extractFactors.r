#' Extracts \code{ArrayData} factors
#' 
#' Extracts the factors, given by an \code{ArrayData} object, that can be used
#' by \code{\link{diffExp}}
#' 
#' 
#' @param arrayData an \code{ArrayData} object.
#' @return A \code{list} with elements:
#' 
#' \item{factors}{Assigns one factor to each array} \item{uniqueFactors}{The
#' unique factors that can be used to form contrasts}
#' @author Leif Varemo \email{piano.rpkg@@gmail.com} and Intawat Nookaew
#' \email{piano.rpkg@@gmail.com}
#' @seealso \pkg{\link{piano}}, \code{\link{diffExp}}
#' @examples
#' 
#'   # Get path to example data and setup files:
#'   dataPath <- system.file("extdata", package="piano")
#' 
#'   # Load normalized data:
#'   myArrayData <- loadMAdata(datadir=dataPath, dataNorm="norm_data.txt.gz", platform="yeast2")
#' 
#'   #Extract the factors that can be used in the call to diffExp:
#'   extractFactors(myArrayData)
#' 
extractFactors <- function(arrayData) {

  # Get factors from setup
  factorAnnot <- arrayData$setup[,1]
  if(dim(arrayData$setup)[2] > 1) {
    for(i in 1:(dim(arrayData$setup)[2]-1)) {
      factorAnnot <- paste(factorAnnot,arrayData$setup[,i+1], sep="_")
    }
  }

  # Order factors according to dataNorm columns
  factors = NA
  for(i in 1:length(colnames(arrayData$dataNorm))) {
    factors[i] <- factorAnnot[rownames(arrayData$setup) == colnames(arrayData$dataNorm)[i]]
  }

  uniqueFactors <- unique(factors)
  factors <- as.data.frame(factors, row.names=colnames(arrayData$dataNorm))
  list(factors=factors,uniqueFactors=uniqueFactors)
}
