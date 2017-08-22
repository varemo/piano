checkLoadArg <- function(statistics,
                         signs, 
                         statMethod, 
                         signMethod, 
                         adjMethod,
                         gsc,
                         gsSize,
                         permStatistics,
                         permSigns,
                         nPerm,
                         gseaParam,
                         ncpus,
                         verbose) {
   
   
   #*********************************************
   #  Check method arguments
   #*********************************************
   
   # Check verbose:
   if(class(verbose) != "logical") stop("argument verbose has to be TRUE or FALSE")
   
   # Check ncpus:
   if(class(ncpus) != "numeric") stop("argument ncpus should be an integer")
   if(length(ncpus) != 1) stop("argument ncpus should be an integer")
   
   # Check statMethod:
   tmp <- try(statMethod <- match.arg(statMethod, c("fisher","stouffer","reporter","tailStrength","wilcoxon",
                                                    "mean","median","sum","maxmean","gsea","fgsea","page"), several.ok=FALSE), silent=TRUE)
   if(class(tmp) == "try-error") {
      stop("argument geneSetStat set to unknown method")
   }
   
   # Check signMethod:
   tmp <- try(signMethod <- match.arg(signMethod, c("geneSampling","samplePermutation","nullDist"), several.ok=FALSE), silent=TRUE)
   if(class(tmp) == "try-error") {
      stop("argument signifMethod set to unknown method")
   }
   if(signMethod == "geneSampling") signMethod <- "geneperm"
   if(signMethod == "samplePermutation") signMethod <- "sampleperm"
   if(signMethod == "nullDist") signMethod <- "distribution"
   
   if(signMethod == "sampleperm" & is.null(permStatistics)) {
      stop("signifMethod='samplePermutation' but no permStats given, can not perform sample permutation")  
   }
   if(signMethod == "geneperm" & !is.null(permStatistics)) {
      warning("signifMethod='geneSampling' and argument permStats given, will not use permStats")  
   }
   if(signMethod == "distribution" & !statMethod %in% c("fisher","stouffer","reporter","wilcoxon","page")) {
      stop(paste("signifMethod='nullDist' is not allowed for geneSetStat='",statMethod,"'",sep=""))
   }
   if(statMethod=="fgsea" & signMethod != "geneperm") {
      stop("only signifMethod='geneSampling' is allowed for geneSetStat='fgsea'") 
   }
   
   # Check adjMethod:
   tmp <- try(adjMethod <- match.arg(adjMethod, c("holm", "hochberg", "hommel", "bonferroni", 
                                                  "BH", "BY","fdr", "none"), several.ok=FALSE), silent=TRUE)
   if(class(tmp) == "try-error") {
      stop("argument adjMethod set to unknown method")
   }
   if(!adjMethod %in% c("fdr","none") & statMethod == "gsea") {
      statMethod <- "fdr"
      warning("adjMethod can only be 'fdr' for statMethod='gsea', using fdr despite different user setting")
   }
   
   # Check nPerm:
   if(length(nPerm) != 1) stop("length of argument nPerm has to be 1")
   if(nPerm < 100) stop("argument nPerm has to be >100")
   if(nPerm%%ncpus != 0) stop("argument ncpus should be set so there is an integer x such that x*ncpus=nPerm")
   
   # Check gseaParam:
   if(gseaParam < 1) stop("gseaParam has to be larger than 0")
   
   
   #*********************************************
   #  Extract statistics and signs
   #*********************************************
   
   # Statistics as matrix like object:
   tmp <- try(statistics <- as.matrix(statistics), silent=TRUE)
   if(class(tmp) == "try-error") {
      stop("argument geneLevelStats could not be converted into a matrix")
   }
   
   # Allow only one comparison, i.e. one-column statistics:
   if(ncol(statistics) != 1) stop("geneLevelStats should be a vector or only contain one column")
   
   # Check for NA:s:
   if(sum(is.na(statistics)) > 0) {
      stop("NA values not allowed in geneLevelStats")  
   }
   
   # Save stat type:
   if(min(statistics) >= 0 & max(statistics) <= 1) {
      statType <- "p"
   } else if(sign(min(statistics)) != sign(max(statistics))) {
      statType <- "t"
   } else {
      statType <- "F"
   }
   
   # Fix zeros and ones in p-values:
   if(statType == "p") {
      statistics[statistics < 1e-100] <- 1e-100
      statistics[statistics > 0.9999999999999999] <- 0.9999999999999999
   }
   
   # Check statistics range:    
   if(statMethod %in% c("fisher","stouffer","reporter","tailStrength") & statType != "p") {
      stop(paste("geneLevelStats does not lie in [0,1], geneLevelStats can only be p-values for geneSetStat='",statMethod,"'",sep=""))
   }
   
   if(statMethod %in% c("maxmean","gsea","fgsea","page") & statType != "t") {
      stop(paste("geneLevelStats has to contain both positive and negative scores for geneSetStat='",statMethod,"'",sep=""))
   }
   
   if(statMethod == "page") {
      for(i in 1:ncol(statistics)) {
         if(stats::sd(statistics[,i]) > 1e300) {
            stop("standard deviation of geneLevelStats is close to infinity, for geneSetStat 'page' all gene-set statistics will be zero, change geneLevelStats or geneSetStat")  
         }
      }
   }
   
   # Check statistics range for other methods:
   # ...no need, wilcoxon, sum, mean, and median can handle all statistic ranges!
   
   # Check for duplicate names in stats (eg 'artifact' of probes => genes):
   # ...now: if duplicate, take both values for gene set.
   tmp <- rownames(statistics)
   if(length(tmp) > length(unique(tmp))) warning("Found duplicates in rownames(geneLevelStats), all values will be used for calculation of gene set statistics. It is recommended to avoid this and handle duplicates prior to running runGSA. In particular, the GSEA and FGSEA implementations will give different results!")
   
   # Count number of genes and conditions in statistics:
   info <- list()
   info$nGenesStatistics <- length(rownames(statistics))
   info$nContrasts <- ncol(statistics)
   
   # Signs as matrix like object:
   if(!is.null(signs) & statType != "t") {
      tmp <- try(signs <- as.matrix(signs), silent=TRUE)
      if(class(tmp) == "try-error") {
         stop("argument directions could not be converted into a matrix")
      }
      
      # Allow only one comparison, i.e. one-column signs:
      if(ncol(signs) != 1) stop("argument directions should be a vector or only contain one column")
      
   # Check if both signs and permSigns are given:
      if(signMethod=="sampleperm" & is.null(permSigns)) {
         stop("for signifMethod='samplePermutation', both directions and permDirections have to be given, or vice versa")  
      }
      
   # Check for zeros: ## Removed this 120615 ##
      #if(sum(signs == 0) != 0) stop("no zeros allowed in argument signs")
   
   # Check for NA:s:   
      if(sum(is.na(signs)) > 0) {
         stop("NA values not allowed in directions")  
      }
      
   # Check signs correlation with statistics:
      if(nrow(statistics) != nrow(signs) | ncol(statistics) != ncol(signs)) {
         stop("dimensions of geneLevelStats and directions are not the same")   
      }
      if(is.null(rownames(statistics)) | is.null(rownames(signs))) {
         stop("rownames of geneLevelStats and directions do not match")  
      }
      if(!identical(rownames(statistics),rownames(signs))) {
         stop("rownames of geneLevelStats and directions do not match")
      }
      
   # Convert to -1/+1:
      signs <- sign(signs)
   
   # Add signs if needed: ## Removed this 120615 ##
      #if(statType %in% c("p","F")) {
      #   statistics <- statistics*signs
      #}
      
   # Update statType:
      if(statType == "p") statType <- "p-signed"
      if(statType == "F") statType <- "F-signed"
      
   } else if(!is.null(signs) & statType == "t") {
      warning("geneLevelStats are t-like and do not require information given by argument directions, argument directions will not be used")  
      signs <- "none"
   } else {
      signs <- "none"  
   }
   
   
   #*********************************************
   #  Data subsets
   #*********************************************
      
#    if(statMethod %in% c("fisher","stouffer","reporter","tailStrength")) {
#       if(statType == "p") {
#          dataSubsets <- c("abs")  
#       } else if(statType == "p-signed" & statMethod != "fisher") {
#          dataSubsets <- c("abs","subset","all")
#       } else { 
#          dataSubsets <- c("abs","subset")      
#       }
#    } else if(statMethod %in% c("mean","median","sum","wilcoxon")) {
#       if(statType %in% c("p","F")) {
#          dataSubsets <- c("abs")  
#       } else if(statType == "p-signed") {
#          dataSubsets <- c("abs","subset","all")
#       } else if(statType == "F-signed") {
#          dataSubsets <- c("abs","subset")
#       } else if(statType == "t") {
#          dataSubsets <- c("all","abs","subset")
#       }
#    } else if(statMethod %in% c("maxmean","gsea","page")) {
#       dataSubsets <- c("all")  
#    }
#    
   
   #*********************************************
   #  Extract gene set collection and additional info
   #*********************************************
   
   # Check if gsc is a list:
   if(class(gsc) != "GSC") stop("argument gsc should have class 'GSC' as output from loadGSC()")

   # Extract and remove addInfo data.frame from gsc:
   addInfo <- gsc$addInfo
      
   # Extract the true gsc part:
   gsc <- gsc$gsc
   
   # Remove genes with no occurance in statistics:
   tmp <- setdiff(unique(unlist(gsc)), rownames(statistics))
   gsc <- lapply(gsc, intersect, rownames(statistics))
   info$removedGenesGSC <- length(tmp)
   
   # Remove gene sets with zero genes:
   tmp <- length(gsc)
   gsc <- gsc[unlist(lapply(gsc,length)) > 0]
   info$removedGSnoGenes <- tmp - length(gsc)
   
   # Remove gene sets not in interval (gsSize[1],gsSize[2]):
   if(length(gsSize) != 2) stop("length(gsSizeLim) has to equal 2")
   if(gsSize[1] > gsSize[2]) stop("gsSizeLim[1] has to be smaller or equal to gsSizeLim[2]")
   if(gsSize[1] < 1) stop("gsSizeLim[1] has to be larger than 0")
   tmp <- length(gsc)
   gsc <- gsc[unlist(lapply(gsc,length)) >= gsSize[1] & unlist(lapply(gsc,length)) <= gsSize[2]]
   info$removedGSsizeLimit <- tmp - length(gsc)
   
   # Number of gene sets with additional info:
   if(class(addInfo) != "character") {
      tmp <- sum(addInfo[,1] %in% names(gsc))
      info$nGeneSetsWithAddInfo <- tmp
   } else {
      info$nGeneSetsWithAddInfo <- 0
   }
   
   # Save info from gsc:
   info$nGenesGSC <- length(unique(unlist(gsc)))
   info$nGeneSets <- length(gsc)

   
   #*********************************************
   #  Extract permStatistics and permSigns
   #*********************************************
   
   # Check that permStatistics is list:
   if(!is.null(permStatistics) & signMethod == "sampleperm") {
      
      permStatistics <- list(permStatistics)
      permSigns <- list(permSigns)
      
      if(!is.null(permSigns) & statType == "t") {
         warning("geneLevelStats are t-like and do not require information given by argument permDirections, argument permDirections will not be used")  
      }
      
      # Check that signs and permSigns co-occur:
      if(is.null(signs) & !is.null(permSigns)) {
         if(statType != "t") stop("for signifMethod='samplePermutation', both directions and permDirections have to be given, or vice versa")  
      }
      
      # For each sample permutation matrix in the list:
      for(i in 1:length(permStatistics)) {
         
         # Sample permutations as matrix like object:
         tmp <- try(permStatistics[[i]] <- as.matrix(permStatistics[[i]]), silent=TRUE)
         if(class(tmp) == "try-error") {
            stop("permStats could not be converted into a matrix")
         }
         
         # Check for NA:s:
         if(sum(is.na(permStatistics[[i]])) > 0) {
            stop("NA values not allowed in permStats")  
         }
         
         # Check stat type:
         if(min(permStatistics[[i]]) >= 0 & max(permStatistics[[i]]) <= 1) {
            tmp <- "p"
         } else if(sign(min(permStatistics[[i]])) != sign(max(permStatistics[[i]]))) {
            tmp <- "t"
         } else {
            tmp <- "F"
         }
         
         # Check same as statistics:
         if(tmp == "p" & !statType%in%c("p","p-signed")) stop("geneLevelStats and permStats has to contain the same type of statistics")
         if(tmp == "F" & !statType%in%c("F","F-signed")) stop("geneLevelStats and permStats has to contain the same type of statistics")
         if(tmp == "t" & !statType == "t") stop("geneLevelStats and permStats has to contain the same type of statistics")
         
         # Fix zeros to small number:
         if(tmp == "p") {
            permStatistics[[i]][permStatistics[[i]] < 1e-100] <- 1e-100
            permStatistics[[i]][permStatistics[[i]] > 0.9999999999999999] <- 0.9999999999999999
         }
         
         # Check statistics range:    
         if(statMethod %in% c("fisher","stouffer","reporter","tailStrength") & tmp != "p") {
            stop(paste("permStats does not lie in [0,1], permStats can only be p-values for geneSetStat='",statMethod,"'",sep=""))
         }
         
         if(statMethod %in% c("maxmean","gsea","page") & tmp != "t") {
            stop(paste("permStats has to contain both positive and negative scores for geneSetStat='",statMethod,"'",sep=""))
         }
         
         # Sample permutations correlation with statistics:
         if(nrow(statistics) != nrow(permStatistics[[i]])) {
            stop("permStats does not have the same number of rows as geneLevelStats")   
         }
         if(!identical(rownames(statistics), rownames(permStatistics[[i]]))) {
            stop("rownames of geneLevelStats and permStats do not match")
         }
         
         # Signs as matrix like object:
         if(!is.null(permSigns) & statType != "t") {
            tmp <- try(permSigns[[i]] <- as.matrix(permSigns[[i]]), silent=TRUE)
            if(class(tmp) == "try-error") {
               stop("permDirections could not be converted into a matrix")
            }
            
            # Check for zeros: ## Removed 120615 ##
            #if(sum(permSigns[[i]] == 0) != 0) stop("no zeros allowed in argument permSigns")
            
            # Check for NA:s:   
            if(sum(is.na(signs[[i]])) > 0) {
               stop("NA values not allowed in permDirections")  
            }
            
            # Check signs correlation with statistics:
            if(nrow(permStatistics[[i]]) != nrow(permSigns[[i]]) | ncol(permStatistics[[i]]) != ncol(permSigns[[i]])) {
               stop("dimensions of permStats and permDirections are not the same")   
            }
            if(!identical(rownames(permStatistics[[i]]),rownames(permSigns[[i]]))) {
               stop("rownames of permStats and permDirections do not match")
            }
            
            # Convert to -1/+1:
            permSigns[[i]] <- sign(permSigns[[i]])
            
            # Add signs if needed: ## Removed 120615 ##
            #permStatistics[[i]] <- permStatistics[[i]]*permSigns[[i]]
            
         } else {
            permSigns <- "none"  
         }
      } 
      
      # Check sample permuation matrix similarity:
      #tmp <- unlist(lapply(permStatistics,ncol))
      #if(!all(tmp == tmp[1])) stop("the number of columns of the elements of list permStatistics differ")
      
      # Check that there exists a matrix for each contrast:
      #if(length(permStatistics) != ncol(statistics)) stop("argument permStatistics has to have the same length as ncol(statistics)")
      
      # Save info:
      info$nSamplePermutations <- ncol(permStatistics[[1]])
      nPerm <- ncol(permStatistics[[1]])
      
   # No permStatistics specified:
   } else {
      permStatistics <- "none"
   }
   
   
   #*********************************************
   #  Return results
   #*********************************************
   
   res <- list()
   res$statistics      <- statistics
   res$statType        <- statType
   res$statMethod      <- statMethod
   res$signMethod      <- signMethod
   res$adjMethod       <- adjMethod
   res$gsc             <- gsc
   res$permStatistics  <- permStatistics
   res$addInfo         <- addInfo
   res$info            <- info
   res$dataSubsets     <- "Not used"
   res$nPerm           <- nPerm
   res$signs           <- signs
   res$permSigns       <- permSigns
   
   return(res)
}
