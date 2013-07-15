runGSA <- function(geneLevelStats,
                   directions=NULL, 
                   geneSetStat="mean", 
                   signifMethod="geneSampling", 
                   adjMethod="fdr",
                   gsc,
                   gsSizeLim=c(1,Inf),
                   permStats=NULL,
                   permDirections=NULL,
                   nPerm=1e4,
                   gseaParam=1,
                   ncpus=1,
                   verbose=TRUE) {
   

   #*********************************************
   #  Check arguments and load data
   #*********************************************

   # Start timing:
   startTime <- proc.time()
   
   if(verbose==TRUE)  {
      message("Running gene set analysis:")
      cat("Checking arguments...")
   }
   tmp <- checkLoadArg(geneLevelStats, directions, geneSetStat, signifMethod, 
                       adjMethod, gsc, gsSizeLim, permStats, permDirections, nPerm, gseaParam, ncpus, verbose)
   
   # Get the output (revert back to original variable names):
   statistics      <- tmp$statistics
   statType        <- tmp$statType 
   statMethod      <- tmp$statMethod 
   signMethod      <- tmp$signMethod 
   adjMethod       <- tmp$adjMethod
   gsc             <- tmp$gsc
   permStatistics  <- tmp$permStatistics
   permSigns       <- tmp$permSigns
   addInfo         <- tmp$addInfo
   info            <- tmp$info
   nPerm           <- tmp$nPerm
   signs           <- tmp$signs
   permSigns       <- tmp$permSigns
   
   if(verbose==TRUE) {
      cat("done!\n")
      if(statMethod=="gsea") message("*** Please note that running the GSEA-method may take a substantial amount of time! ***")
      if(statMethod=="wilcoxon") message("*** Please note that running the Wilcoxon-method may take a substantial amount of time! ***")
      message(paste("Final gene/gene-set association:",info$nGenesGSC,"genes and",info$nGeneSets,"gene sets"))
      message("  Details:")
      message(paste("  Calculating gene set statistics from ",info$nGenesGSC," out of ",info$nGenesStatistics," gene-level statistics",sep=""))
      if(signMethod %in% c("geneperm","sampleperm")) message(paste("  Using all ",info$nGenesStatistics," gene-level statistics for significance estimation",sep=""))
      #message(paste("  Removed",info$redundantGSCrows,"redundant gene/gene-set associations"))
      message(paste("  Removed",info$removedGenesGSC,"genes from GSC due to lack of matching gene statistics"))
      message(paste("  Removed",info$removedGSnoGenes,"gene sets containing no genes after gene removal"))
      message(paste("  Removed additionally",info$removedGSsizeLimit,"gene sets not matching the size limits"))
      message(paste("  Loaded additional information for",info$nGeneSetsWithAddInfo,"gene sets"))
   }   
   
   
   #*********************************************
   #  Calc gs stat
   #*********************************************
   
   if(verbose==TRUE) cat("Calculating gene set statistics...")
   tmp <- GSCstatBatch(statistics, statType, gsc, statMethod, signMethod, gseaParam, signs)
   
   # Get number of genes in each gene set:
   nGenes   <- tmp$nGenes
   nGenesUp <- tmp$nGenesUp
   nGenesDn <- tmp$nGenesDn
   
   # Get the gene set statistics:
   gsStatsAll <- tmp$statsAll
   gsStatsAllTestUp <- tmp$statsAllTestUp
   gsStatsAllTestDn <- tmp$statsAllTestDn
   gsStatsAbs <- tmp$statsAbs
   gsStatsUp  <- tmp$statsUp
   gsStatsDn  <- tmp$statsDn
   
   # Get gene-set stat name:
   gsStatName <- tmp$statName
   
   if(verbose==TRUE) cat("done!\n")
   
   #*********************************************
   #  Calc gs significance
   #*********************************************
   
   if(verbose==TRUE) cat("Calculating gene set significance...")
   
   if(!(statMethod == "wilcoxon" & signMethod == "distribution")) {
   tmp <- GSCsignificanceBatch(statistics, statType, signs, gsc, statMethod, signMethod, permStatistics, permSigns,
                               nGenes, nGenesUp, nGenesDn, gsStatsAll, gsStatsAllTestUp, gsStatsAllTestDn, gsStatsAbs, 
                               gsStatsUp, gsStatsDn, nPerm, gseaParam, ncpus)
   
   }
   
   if(verbose==TRUE) cat("done!\n")
   
   # Get the p-values:
   pValuesAll   <- tmp$pValuesAll
   pValuesAllUp <- tmp$pValuesAllUp
   pValuesAllDn <- tmp$pValuesAllDn
   pValuesAbs   <- tmp$pValuesAbs
   pValuesUp    <- tmp$pValuesUp
   pValuesDn    <- tmp$pValuesDn
   
   # For GSEA:
   gsStatsAllPerm <- tmp$gsStatsAllPerm
   
   
   #*********************************************
   #  Adj. for multiple testing
   #*********************************************
   
   if(verbose==TRUE & adjMethod != "none") cat("Adjusting for multiple testing...")
   
   if(statMethod == "gsea" & adjMethod != "none") {
      tmp <- fdrGSEA(gsStatsAll,gsStatsAllPerm,nGenes,signMethod)
      pValuesAllUpAdj <- tmp$pValuesAllUpAdj
      pValuesAllDnAdj <- tmp$pValuesAllDnAdj
      pValuesAllAdj <- pValuesAllUpAdj*NA
      pValuesAbsAdj <- pValuesAllUpAdj*NA
      pValuesUpAdj <- pValuesAllUpAdj*NA
      pValuesDnAdj <- pValuesAllUpAdj*NA
   } else {
      pValuesAllAdj <- apply(pValuesAll,2,p.adjust,method=adjMethod)
      pValuesAllUpAdj <- apply(pValuesAllUp,2,p.adjust,method=adjMethod)
      pValuesAllDnAdj <- apply(pValuesAllDn,2,p.adjust,method=adjMethod)
      pValuesAbsAdj <- apply(pValuesAbs,2,p.adjust,method=adjMethod)
      pValuesUpAdj <- apply(pValuesUp,2,p.adjust,method=adjMethod)
      pValuesDnAdj <- apply(pValuesDn,2,p.adjust,method=adjMethod)
   }
   
   if(verbose==TRUE & adjMethod != "none") cat("done!\n")
   
   
   #*********************************************
   #  Return results
   #*********************************************
   
   res <- list()
   
   # General info:
   res$geneStatType  <- statType 
   res$geneSetStat    <- statMethod 
   if(signMethod == "geneperm") res$signifMethod <- "geneSampling"
   if(signMethod == "sampleperm") res$signifMethod <- "samplePermutation"
   if(signMethod == "distribution") res$signifMethod <- "nullDist"
   res$adjMethod     <- adjMethod
   res$info          <- info
   #res$dataSubsets  <- dataSubsets
   res$gsSizeLim        <- gsSizeLim
   res$gsStatName    <- gsStatName
   res$nPerm         <- nPerm
   res$gseaParam     <- gseaParam
   #res$contrastName  <- colnames(statistics)
   
   # Gene-level statistics and signs:
   res$geneLevelStats <- statistics
   res$directions <- signs
   
   # GSC info:
   res$gsc          <- gsc
   res$nGenesTot       <- nGenes
   res$nGenesUp  <- nGenesUp
   res$nGenesDn  <- nGenesDn
   
   # Gene set statistics:
   colnames(gsStatsAll)       <- colnames(statistics)
   colnames(gsStatsAllTestUp) <- colnames(statistics)
   colnames(gsStatsAllTestDn) <- colnames(statistics)
   colnames(gsStatsAbs)       <- colnames(statistics)
   colnames(gsStatsUp)        <- colnames(statistics)
   colnames(gsStatsDn)        <- colnames(statistics)
   
   res$statDistinctDir   <- signif(gsStatsAll,5)
   res$statDistinctDirUp <- signif(gsStatsAllTestUp,5)
   res$statDistinctDirDn <- signif(gsStatsAllTestDn,5)
   res$statNonDirectional   <- signif(gsStatsAbs,5)
   res$statMixedDirUp <- signif(gsStatsUp,5)
   res$statMixedDirDn <- signif(gsStatsDn,5)
   
   # Gene set p-values:
   #colnames(pValuesAll)   <- colnames(statistics)
   colnames(pValuesAllUp) <- colnames(statistics)
   colnames(pValuesAllDn) <- colnames(statistics)
   colnames(pValuesAbs)   <- colnames(statistics)
   colnames(pValuesUp)    <- colnames(statistics)
   colnames(pValuesDn)    <- colnames(statistics)
   
   #res$pValuesAll   <- pValuesAll
   res$pDistinctDirUp <- signif(pValuesAllUp,5)
   res$pDistinctDirDn <- signif(pValuesAllDn,5)
   res$pNonDirectional   <- signif(pValuesAbs,5)
   res$pMixedDirUp <- signif(pValuesUp,5)
   res$pMixedDirDn <- signif(pValuesDn,5)
   
   # Gene set adjusted p-values:
   #colnames(pValuesAllAdj)   <- colnames(statistics)  
   colnames(pValuesAllUpAdj) <- colnames(statistics)
   colnames(pValuesAllDnAdj) <- colnames(statistics)
   colnames(pValuesAbsAdj)   <- colnames(statistics)
   colnames(pValuesUpAdj)    <- colnames(statistics)
   colnames(pValuesDnAdj)    <- colnames(statistics)
   
   #res$pValuesAllAdj   <- pValuesAllAdj
   res$pAdjDistinctDirUp <- signif(pValuesAllUpAdj,5)
   res$pAdjDistinctDirDn <- signif(pValuesAllDnAdj,5)
   res$pAdjNonDirectional  <- signif(pValuesAbsAdj,5)
   res$pAdjMixedDirUp <- signif(pValuesUpAdj,5)
   res$pAdjMixedDirDn <- signif(pValuesDnAdj,5)
   
   # Add time info:
   res$runtime <- proc.time() - startTime
   
   class(res) <- "GSAres"
   return(res)
   
}