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