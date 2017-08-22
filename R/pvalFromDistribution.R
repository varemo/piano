pvalFromDistribution <- function(statType,statMethod,nGenes,nGenesUp,nGenesDn,gsStatsAll,gsStatsAllTestUp,gsStatsAllTestDn,gsStatsAbs,gsStatsUp,gsStatsDn,
                                 gsStatsAllPerm,gsStatsAllTestUpPerm,gsStatsAllTestDnPerm,gsStatsAbsPerm,gsStatsUpPerm,gsStatsDnPerm) {
   
   
   # Create result list object:
   res <- list()
   res$pValuesAll <- vector()
   res$pValuesAllUp <- vector()
   res$pValuesAllDn <- vector()
   res$pValuesAbs <- vector()
   res$pValuesUp <- vector()
   res$pValuesDn <- vector()
   
   # Loop over contrasts:
   for(iContrast in 1:ncol(nGenes)) {
      
      # Preallocate:
      pValuesAll <- rep(NA,nrow(nGenes))
      pValuesAllUp <- rep(NA,nrow(nGenes))
      pValuesAllDn <- rep(NA,nrow(nGenes))
      pValuesAbs <- rep(NA,nrow(nGenes))
      pValuesUp <- rep(NA,nrow(nGenes))
      pValuesDn <- rep(NA,nrow(nGenes))
      
      for(iGeneSet in 1:nrow(nGenes)) {
         
         #************************  
         # Fisher:
         #************************
         if(statMethod == "fisher") {
            
            # "Mixed":
            stat <- gsStatsAbs[iGeneSet,iContrast]
            param <- 2*nGenes[iGeneSet,iContrast]
            pValuesAbs[iGeneSet] <- pchisq(stat, param, lower.tail = FALSE)
            
            # Subset up:
            if(statType == "p-signed") {
               stat <- gsStatsUp[iGeneSet,iContrast]
               param <- 2*nGenesUp[iGeneSet,iContrast]
               pValuesUp[iGeneSet] <- pchisq(stat, param, lower.tail = FALSE)
               
               # Subset dn:
               stat <- gsStatsDn[iGeneSet,iContrast]
               param <- 2*nGenesDn[iGeneSet,iContrast]
               pValuesDn[iGeneSet] <- pchisq(stat, param, lower.tail = FALSE)
            }
         }
         
         #************************  
         # Stouffer:
         #************************
         if(statMethod == "stouffer") {
            
            # "Directional":
            if(statType == "p-signed") {
               stat <- gsStatsAllTestUp[iGeneSet,iContrast]
               pValuesAllUp[iGeneSet] <- pnorm(stat,lower.tail=FALSE)
               stat <- gsStatsAllTestDn[iGeneSet,iContrast]
               pValuesAllDn[iGeneSet] <- pnorm(stat,lower.tail=FALSE)
            }
            
            # "Mixed":
            stat <- gsStatsAbs[iGeneSet,iContrast]
            pValuesAbs[iGeneSet] <- pnorm(stat,lower.tail=FALSE)
            
            # Subset up:
            if(statType == "p-signed") {
               stat <- gsStatsUp[iGeneSet,iContrast]
               pValuesUp[iGeneSet] <- pnorm(stat,lower.tail=FALSE)
               
               # Subset dn:
               stat <- gsStatsDn[iGeneSet,iContrast]
               pValuesDn[iGeneSet] <- pnorm(stat,lower.tail=FALSE)
            }
         }
         
         #************************  
         # Reporter:
         #************************
         if(statMethod == "reporter") {
            
            # "Directional":
            if(statType == "p-signed") {
               nIndex <- nGenes[iGeneSet,iContrast]
               
               backgroundTestUp <- gsStatsAllTestUpPerm[[iContrast]][as.character(nIndex),]
               zNormTestUp <- (gsStatsAllTestUp[iGeneSet,iContrast] - mean(backgroundTestUp))/stats::sd(backgroundTestUp)
               pValuesAllUp[iGeneSet] <- 1 - pnorm(zNormTestUp)
               
               backgroundTestDn <- gsStatsAllTestDnPerm[[iContrast]][as.character(nIndex),]
               zNormTestDn <- (gsStatsAllTestDn[iGeneSet,iContrast] - mean(backgroundTestDn))/stats::sd(backgroundTestDn)
               pValuesAllDn[iGeneSet] <- 1 - pnorm(zNormTestDn)
            }
            
            # Mixed:
            nIndex <- nGenes[iGeneSet,iContrast]
            background <- gsStatsAbsPerm[[iContrast]][as.character(nIndex),]
            zNorm <- (gsStatsAbs[iGeneSet,iContrast] - mean(background))/stats::sd(background)
            pValuesAbs[iGeneSet] <- 1 - pnorm(zNorm)
            
            # Subset up:
            if(statType == "p-signed") {
               nIndex <- nGenesUp[iGeneSet,iContrast]
               background <- gsStatsUpPerm[[iContrast]][as.character(nIndex),]
               zNorm <- (gsStatsUp[iGeneSet,iContrast] - mean(background))/stats::sd(background)
               pValuesUp[iGeneSet] <- 1 - pnorm(zNorm)
               
               # Subset dn:
               nIndex <- nGenesDn[iGeneSet,iContrast]
               background <- gsStatsDnPerm[[iContrast]][as.character(nIndex),]
               zNorm <- (gsStatsDn[iGeneSet,iContrast] - mean(background))/stats::sd(background)
               pValuesDn[iGeneSet] <- 1 - pnorm(zNorm)
            }
         }
         
         #************************  
         # Page:
         #************************
         if(statMethod == "page") {
            
            # Directional:
            stat <- gsStatsAll[iGeneSet,iContrast]
            pValuesAllUp[iGeneSet] <- 1 - pnorm(stat)
            pValuesAllDn[iGeneSet] <- 1 - pValuesAllUp[iGeneSet]
            # Old version: pValuesAll[iGeneSet] <- pnorm(stat,lower.tail=stat<0)*2
            # PGSEA: p.val <- 2 * pnorm(-abs(stat2))
         }
      }
      
      # Save results:
      res$pValuesAll <- cbind(res$pValuesAll, pValuesAll)
      res$pValuesAllUp <- cbind(res$pValuesAllUp, pValuesAllUp)
      res$pValuesAllDn <- cbind(res$pValuesAllDn, pValuesAllDn)
      res$pValuesAbs <- cbind(res$pValuesAbs, pValuesAbs)
      res$pValuesUp <- cbind(res$pValuesUp, pValuesUp)
      res$pValuesDn <- cbind(res$pValuesDn, pValuesDn)
   }
   
   return(res)
}