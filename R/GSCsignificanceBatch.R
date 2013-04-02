GSCsignificanceBatch <- function(statistics, statType, signs, gsc, statMethod, signMethod, permStatistics, permSigns,
                                 nGenes, nGenesUp, nGenesDn, gsStatsAll, gsStatsAllTestUp, gsStatsAllTestDn, gsStatsAbs, gsStatsUp, gsStatsDn, nPerm, gseaParam) {
   
   #*********************************************
   # Gene permutation:
   #*********************************************
   
   if(signMethod == "geneperm" | (statMethod == "reporter" & signMethod == "distribution") ) {
      
      # Calculate gene set statistics distributions for each size:
      res <- GSCstatGenePerm(statistics, signs, gsc, statType, statMethod, nGenes, nGenesUp, nGenesDn, nPerm, gseaParam)
      gsStatsAllPerm <- res$gsStatsAllPerm
      gsStatsAllTestUpPerm <- res$gsStatsAllTestUpPerm
      gsStatsAllTestDnPerm <- res$gsStatsAllTestDnPerm
      gsStatsAbsPerm <- res$gsStatsAbsPerm
      gsStatsUpPerm <- res$gsStatsUpPerm
      gsStatsDnPerm <- res$gsStatsDnPerm
      
      
   #*********************************************
   # Sample permutation:   
   #*********************************************
   
   } else if(signMethod == "sampleperm") {
      res <- GSCstatSamplePerm(permStatistics, permSigns, statType, gsc, statMethod, gseaParam)
      gsStatsAllPerm <- res$gsStatsAllPerm
      gsStatsAllTestUpPerm <- res$gsStatsAllTestUpPerm
      gsStatsAllTestDnPerm <- res$gsStatsAllTestDnPerm
      gsStatsAbsPerm <- res$gsStatsAbsPerm
      gsStatsUpPerm <- res$gsStatsUpPerm
      gsStatsDnPerm <- res$gsStatsDnPerm
   
   } else {
      gsStatsAllPerm <- NA
      gsStatsAllTestUpPerm <- NA
      gsStatsAllTestDnPerm <- NA
      gsStatsAbsPerm <- NA
      gsStatsUpPerm <- NA
      gsStatsDnPerm <- NA
   }
   
   
   #*********************************************
   # p-value from fraction
   #*********************************************
      
   if(signMethod == "geneperm") {
      res <- pvalFromFractionGenePerm(statType,statMethod,nGenes,nGenesUp,nGenesDn,gsStatsAll,gsStatsAllTestUp,gsStatsAllTestDn,gsStatsAbs,gsStatsUp,gsStatsDn,
                              gsStatsAllPerm,gsStatsAllTestUpPerm,gsStatsAllTestDnPerm,gsStatsAbsPerm,gsStatsUpPerm,gsStatsDnPerm)
   } else if(signMethod == "sampleperm") {
      res <- pvalFromFractionSamplePerm(statType,statMethod,nGenes,nGenesUp,nGenesDn,gsStatsAll,gsStatsAllTestUp,gsStatsAllTestDn,gsStatsAbs,gsStatsUp,gsStatsDn,
                              gsStatsAllPerm,gsStatsAllTestUpPerm,gsStatsAllTestDnPerm,gsStatsAbsPerm,gsStatsUpPerm,gsStatsDnPerm)
   }
   
   #*********************************************
   # p-value from distribution
   #*********************************************
   
   if(signMethod == "distribution") {
      res <- pvalFromDistribution(statType,statMethod,nGenes,nGenesUp,nGenesDn,gsStatsAll,gsStatsAllTestUp,gsStatsAllTestDn,gsStatsAbs,gsStatsUp,gsStatsDn,
                                  gsStatsAllPerm,gsStatsAllTestUpPerm,gsStatsAllTestDnPerm,gsStatsAbsPerm,gsStatsUpPerm,gsStatsDnPerm)
   }
   
   #*********************************************
   # Save results:
   #*********************************************
   res$gsStatsAllPerm <- gsStatsAllPerm # For GSEA FDR calculation
   return(res)
}