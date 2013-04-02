pvalFromFractionGenePerm <- function(statType,statMethod,nGenes,nGenesUp,nGenesDn,gsStatsAll,gsStatsAllTestUp,gsStatsAllTestDn,gsStatsAbs,gsStatsUp,gsStatsDn,
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
            value <- gsStatsAbs[iGeneSet,iContrast]  
            nIndex <- nGenes[iGeneSet,iContrast]
            background <- gsStatsAbsPerm[[iContrast]][as.character(nIndex),]
            pValuesAbs[iGeneSet] <- sum(background >= value)/length(background)
            
            # Subset up:
            if(statType == "p-signed") {
               nIndex <- nGenesUp[iGeneSet,iContrast]
               if(nIndex != 0) {
                  value <- gsStatsUp[iGeneSet,iContrast]  
                  background <- gsStatsUpPerm[[iContrast]][as.character(nIndex),]
                  pValuesUp[iGeneSet] <- sum(background >= value)/length(background)
               }
               
               # Subset dn:
               nIndex <- nGenesDn[iGeneSet,iContrast]
               if(nIndex != 0) {
                  value <- gsStatsDn[iGeneSet,iContrast]  
                  background <- gsStatsDnPerm[[iContrast]][as.character(nIndex),]
                  pValuesDn[iGeneSet] <- sum(background >= value)/length(background)
               }
            }
         }
         
         
         #************************
         # Stouffer:
         # Reporter:
         # Tail strength:
         #************************
         if(statMethod %in% c("stouffer","reporter","tailStrength")) {

            # "Directional":
            if(statType == "p-signed") {
               valueTestUp <- gsStatsAllTestUp[iGeneSet,iContrast]
               valueTestDn <- gsStatsAllTestDn[iGeneSet,iContrast]
               nIndex <- nGenes[iGeneSet,iContrast]   
               backgroundTestUp <- gsStatsAllTestUpPerm[[iContrast]][as.character(nIndex),]
               backgroundTestDn <- gsStatsAllTestDnPerm[[iContrast]][as.character(nIndex),]
               pValuesAllUp[iGeneSet] <- sum(backgroundTestUp >= valueTestUp)/length(backgroundTestUp) # all up
               pValuesAllDn[iGeneSet] <- sum(backgroundTestDn >= valueTestDn)/length(backgroundTestDn) # all dn
            }
            
            # "Mixed":
            value <- gsStatsAbs[iGeneSet,iContrast]  
            nIndex <- nGenes[iGeneSet,iContrast]
            background <- gsStatsAbsPerm[[iContrast]][as.character(nIndex),]
            pValuesAbs[iGeneSet] <- sum(background >= value)/length(background)
            
            # Subset up:
            if(statType == "p-signed") {
               nIndex <- nGenesUp[iGeneSet,iContrast]
               if(nIndex != 0) {
                  value <- gsStatsUp[iGeneSet,iContrast]  
                  background <- gsStatsUpPerm[[iContrast]][as.character(nIndex),]
                  pValuesUp[iGeneSet] <- sum(background >= value)/length(background)
               }
               
               # Subset dn:
               nIndex <- nGenesDn[iGeneSet,iContrast]
               if(nIndex != 0) {
                  value <- gsStatsDn[iGeneSet,iContrast]  
                  background <- gsStatsDnPerm[[iContrast]][as.character(nIndex),]
                  pValuesDn[iGeneSet] <- sum(background >= value)/length(background)
               }
            }
         }
         
         
         #************************  
         # Wilcoxon:
         #************************
         if(statMethod == "wilcoxon") {
            
            # "Directional":
            nIndex <- nGenes[iGeneSet,iContrast]
            if(statType == "t") {
               value <- gsStatsAll[iGeneSet,iContrast]  
               background <- gsStatsAllPerm[[iContrast]][as.character(nIndex),]
               pValuesAllUp[iGeneSet] <- sum(background >= value)/length(background) # all up
               pValuesAllDn[iGeneSet] <- sum(background <= value)/length(background) # all dn
            } else if(statType == "p-signed") {
               valueTestUp <- gsStatsAllTestUp[iGeneSet,iContrast]
               valueTestDn <- gsStatsAllTestDn[iGeneSet,iContrast]
               backgroundTestUp <- gsStatsAllTestUpPerm[[iContrast]][as.character(nIndex),]
               backgroundTestDn <- gsStatsAllTestDnPerm[[iContrast]][as.character(nIndex),]
               pValuesAllUp[iGeneSet] <- sum(backgroundTestUp <= valueTestUp)/length(backgroundTestUp) # all up
               pValuesAllDn[iGeneSet] <- sum(backgroundTestDn <= valueTestDn)/length(backgroundTestDn) # all dn
            }
            
            # "Mixed":
            value <- gsStatsAbs[iGeneSet,iContrast]  
            nIndex <- nGenes[iGeneSet,iContrast]
            background <- gsStatsAbsPerm[[iContrast]][as.character(nIndex),]
            if(statType %in% c("p","p-signed")) {
               pValuesAbs[iGeneSet] <- sum(background <= value)/length(background)
            } else {
               pValuesAbs[iGeneSet] <- sum(background >= value)/length(background)
            }
            
            # Subset up:
            nIndex <- nGenesUp[iGeneSet,iContrast]
            if(nIndex != 0) {
               value <- gsStatsUp[iGeneSet,iContrast]  
               background <- gsStatsUpPerm[[iContrast]][as.character(nIndex),]
               if(statType == "p-signed") {
                  pValuesUp[iGeneSet] <- sum(background <= value)/length(background)
               } else if(statType %in% c("t","F-signed")) {
                  pValuesUp[iGeneSet] <- sum(background >= value)/length(background)
               }
            }
            
            # Subset dn:
            nIndex <- nGenesDn[iGeneSet,iContrast]
            if(nIndex != 0) {
               value <- gsStatsDn[iGeneSet,iContrast]  
               background <- gsStatsDnPerm[[iContrast]][as.character(nIndex),]
               if(statType == "p-signed") {
                  pValuesDn[iGeneSet] <- sum(background <= value)/length(background)
               } else if(statType %in% c("t","F-signed")) {
                  pValuesDn[iGeneSet] <- sum(background >= value)/length(background)
               }
            }  
         }
         
         #************************
         # Mean:
         # Median:
         # Sum:
         #************************
         if(statMethod %in% c("mean","median","sum")) {
            
            # "Directional":
            nIndex <- nGenes[iGeneSet,iContrast]
            if(statType == "t") {
               value <- gsStatsAll[iGeneSet,iContrast]
               background <- gsStatsAllPerm[[iContrast]][as.character(nIndex),]
               pValuesAllUp[iGeneSet] <- sum(background >= value)/length(background) # all up
               pValuesAllDn[iGeneSet] <- sum(background <= value)/length(background) # all dn
            } else if(statType == "p-signed") {
               valueTestUp <- gsStatsAllTestUp[iGeneSet,iContrast]
               valueTestDn <- gsStatsAllTestDn[iGeneSet,iContrast]
               backgroundTestUp <- gsStatsAllTestUpPerm[[iContrast]][as.character(nIndex),]
               backgroundTestDn <- gsStatsAllTestDnPerm[[iContrast]][as.character(nIndex),]
               pValuesAllUp[iGeneSet] <- sum(backgroundTestUp <= valueTestUp)/length(backgroundTestUp) # all up
               pValuesAllDn[iGeneSet] <- sum(backgroundTestDn <= valueTestDn)/length(backgroundTestDn) # all dn
            }
            
            # "Mixed":
            value <- gsStatsAbs[iGeneSet,iContrast]  
            nIndex <- nGenes[iGeneSet,iContrast]
            background <- gsStatsAbsPerm[[iContrast]][as.character(nIndex),]
            if(statType %in% c("p","p-signed")) {
               pValuesAbs[iGeneSet] <- sum(background <= value)/length(background)
            } else {
               pValuesAbs[iGeneSet] <- sum(background >= value)/length(background)
            }
            
            # Subset up:
            nIndex <- nGenesUp[iGeneSet,iContrast]
            if(nIndex != 0) {
               value <- gsStatsUp[iGeneSet,iContrast]  
               background <- gsStatsUpPerm[[iContrast]][as.character(nIndex),]
               if(statType == "p-signed") {
                  pValuesUp[iGeneSet] <- sum(background <= value)/length(background)
               } else if(statType %in% c("t","F-signed")){
                  pValuesUp[iGeneSet] <- sum(background >= value)/length(background)
               }
            }
            
            # Subset dn:
            nIndex <- nGenesDn[iGeneSet,iContrast]
            if(nIndex != 0) {
               value <- gsStatsDn[iGeneSet,iContrast]  
               background <- gsStatsDnPerm[[iContrast]][as.character(nIndex),]
               if(statType == "p-signed") {
                  pValuesDn[iGeneSet] <- sum(background <= value)/length(background)
               } else if(statType %in% c("t","F-signed")){
                  pValuesDn[iGeneSet] <- sum(background >= value)/length(background)
               }
            }
         }
         
         #************************   
         # Maxmean:
         #************************
         if(statMethod == "maxmean") {
            
            # "Mixed":
            value <- gsStatsAbs[iGeneSet,iContrast]  
            nIndex <- nGenes[iGeneSet,iContrast]
            background <- gsStatsAbsPerm[[iContrast]][as.character(nIndex),]
            pValuesAbs[iGeneSet] <- sum(background >= value)/length(background)  
         }
         
         #************************   
         # GSEA:
         #************************
         if(statMethod == "gsea") {
            value <- gsStatsAll[iGeneSet,iContrast]  
            nIndex <- nGenes[iGeneSet,iContrast]   
            background <- gsStatsAllPerm[[iContrast]][as.character(nIndex),]
            if(value < 0) {
               backgroundNeg <- background[background < 0]
               pValuesAllDn[iGeneSet] <- sum(backgroundNeg <= value)/length(backgroundNeg) # all
               pValuesAllUp[iGeneSet] <- NA
            } else {
               backgroundPos <- background[background > 0]
               pValuesAllUp[iGeneSet] <- sum(backgroundPos >= value)/length(backgroundPos) # all
               pValuesAllDn[iGeneSet] <- NA
            }
            
         }
         
         #************************  
         # Page:
         #************************
         if(statMethod == "page") {
            value <- gsStatsAll[iGeneSet,iContrast]  
            nIndex <- nGenes[iGeneSet,iContrast]   
            background <- gsStatsAllPerm[[iContrast]][as.character(nIndex),]
            pValuesAllUp[iGeneSet] <- sum(background >= value)/length(background) # all up
            pValuesAllDn[iGeneSet] <- sum(background <= value)/length(background) # all dn
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