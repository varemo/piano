print.ArrayData <- function(x, ...){

  if((length(x) == 4 | length(x) == 3 | length(x) == 2) & "dataNorm" %in% attributes(x)$names & "setup" %in% attributes(x)$names) {
    cat("ArrayData object\n")

    if("dataRaw" %in% attributes(x)$names) {
      cat("  dataRaw:    raw data as an AffyBatch object\n")
      cat("              samples: ")
      cat(length(sampleNames(x$dataRaw)))
      cat("\n")
      cat("              genes: ")
      cat(length(featureNames(x$dataRaw)))
      cat("\n")
    }
    cat("  dataNorm:   data frame containing normalized expression values\n")
    cat("              samples: ")
    cat(dim(x$dataNorm)[2])
    cat("\n")
    cat("              genes: ")
    cat(dim(x$dataNorm)[1])
    cat("\n")
    cat("  setup:      data frame containing experimental setup\n")
    cat("              samples: ")
    cat(dim(x$setup)[1])
    cat("\n")
    cat("              factor categories: ")
    cat(dim(x$setup)[2])
    cat("\n")
    if("annotation" %in% attributes(x)$names) {
      cat("  annotation: data frame containing annotation information\n")
      cat("              genes: ")
      cat(dim(x$annotation)[1])
      cat("\n")
    }
 
  } else {
    stop("this object contains errors")
  }
  
}