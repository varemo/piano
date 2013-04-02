GSCstatSamplePerm <- function(permStatistics, permSigns, statType, gsc, method, gseaParam) {
   
   # Create result list object:
   res <- list()
   res$gsStatsAllPerm <- list()
   res$gsStatsAllTestUpPerm <- list()
   res$gsStatsAllTestDnPerm <- list()
   res$gsStatsAbsPerm <- list()
   res$gsStatsUpPerm <- list()
   res$gsStatsDnPerm <- list()
   for(iContrast in 1:length(permStatistics)) {
      res$gsStatsAllPerm[[iContrast]] <- vector()
      res$gsStatsAllTestUpPerm[[iContrast]] <- vector()
      res$gsStatsAllTestDnPerm[[iContrast]] <- vector()
      res$gsStatsAbsPerm[[iContrast]] <- vector()
      res$gsStatsUpPerm[[iContrast]] <- vector()
      res$gsStatsDnPerm[[iContrast]] <- vector()
   }
   
   # For each sample permutation, calculate new gene-set statistics for all contrasts:
   for(iPerm in 1:ncol(permStatistics[[1]])) {
      
      # Get statistics and signs for all contrast for permutation iPerm:
      statistics <- vector()
      signs <- vector()
      for(iContrast in 1:length(permStatistics)) {
         statistics <- cbind(statistics,permStatistics[[iContrast]][,iPerm])
         if(statType %in% c("p-signed","F-signed")) {
            signs <- cbind(signs,permSigns[[iContrast]][,iPerm])
         } else {
            signs <- "none"
         }
      }
      
      # Calculate gene-set statistics (for permutation iPerm):
      if(method == "wilcoxon") method <- "wilcoxon_fast"
      tmp <- GSCstatBatch(statistics, statType, gsc, method, "sampleperm", gseaParam, signs) # <------ fast wilcoxon, added!
      gsStatsAll <- tmp$statsAll
      gsStatsAllTestUp <- tmp$statsAllTestUp
      gsStatsAllTestDn <- tmp$statsAllTestDn
      gsStatsAbs <- tmp$statsAbs
      gsStatsUp  <- tmp$statsUp
      gsStatsDn  <- tmp$statsDn
      
      # Save results:
      for(iContrast in 1:length(permStatistics)) {
         res$gsStatsAllPerm[[iContrast]] <- cbind(res$gsStatsAllPerm[[iContrast]],gsStatsAll[,iContrast])
         res$gsStatsAllTestUpPerm[[iContrast]] <- cbind(res$gsStatsAllTestUpPerm[[iContrast]],gsStatsAllTestUp[,iContrast])
         res$gsStatsAllTestDnPerm[[iContrast]] <- cbind(res$gsStatsAllTestDnPerm[[iContrast]],gsStatsAllTestDn[,iContrast])
         res$gsStatsAbsPerm[[iContrast]] <- cbind(res$gsStatsAbsPerm[[iContrast]],gsStatsAbs[,iContrast])
         res$gsStatsUpPerm[[iContrast]] <- cbind(res$gsStatsUpPerm[[iContrast]],gsStatsUp[,iContrast])
         res$gsStatsDnPerm[[iContrast]] <- cbind(res$gsStatsDnPerm[[iContrast]],gsStatsDn[,iContrast])
      }
   }
   
   return(res)
   
}