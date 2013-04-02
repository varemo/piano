pvalFromFractionSamplePerm <- function(statType,statMethod,nGenes,nGenesUp,nGenesDn,gsStatsAll,gsStatsAllTestUp,gsStatsAllTestDn,gsStatsAbs,gsStatsUp,gsStatsDn,
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
            background <- gsStatsAbsPerm[[iContrast]][iGeneSet,]
            pValuesAbs[iGeneSet] <- sum(background >= value, na.rm=TRUE)/sum(!is.na(background))
            
#             # Subset up:
#             if(statType == "p-signed") {
#                if(!is.na(nGenesUp[iGeneSet,iContrast])) {
#                   if(nGenesUp[iGeneSet,iContrast] != 0) {
#                      value <- gsStatsUp[iGeneSet,iContrast]  
#                      background <- gsStatsUpPerm[[iContrast]][iGeneSet,]
#                      pValuesUp[iGeneSet] <- sum(background >= value, na.rm=TRUE)/sum(!is.na(background))
#                   }
#                }
#                
#                # Subset dn:
#                if(!is.na(nGenesDn[iGeneSet,iContrast])) {
#                   if(nGenesDn[iGeneSet,iContrast] != 0) {
#                      value <- gsStatsDn[iGeneSet,iContrast]  
#                      background <- gsStatsDnPerm[[iContrast]][iGeneSet,]
#                      pValuesDn[iGeneSet] <- sum(background >= value, na.rm=TRUE)/sum(!is.na(background))
#                   }
#                }
#             }
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
                 
               backgroundTestUp <- gsStatsAllTestUpPerm[[iContrast]][iGeneSet,]
               backgroundTestDn <- gsStatsAllTestDnPerm[[iContrast]][iGeneSet,]
               
               pValuesAllUp[iGeneSet] <- sum(backgroundTestUp >= valueTestUp, na.rm=TRUE)/sum(!is.na(backgroundTestUp)) # all up
               pValuesAllDn[iGeneSet] <- sum(backgroundTestDn >= valueTestDn, na.rm=TRUE)/sum(!is.na(backgroundTestDn)) # all dn
            }
            
            # "Mixed":
            value <- gsStatsAbs[iGeneSet,iContrast]  
            background <- gsStatsAbsPerm[[iContrast]][iGeneSet,]
            pValuesAbs[iGeneSet] <- sum(background >= value, na.rm=TRUE)/sum(!is.na(background))
            
#             # Subset up:
#             if(statType == "p-signed") {
#                if(!is.na(nGenesUp[iGeneSet,iContrast])) {
#                   if(nGenesUp[iGeneSet,iContrast] != 0) {
#                      value <- gsStatsUp[iGeneSet,iContrast]  
#                      background <- gsStatsUpPerm[[iContrast]][iGeneSet,]
#                      pValuesUp[iGeneSet] <- sum(background >= value, na.rm=TRUE)/sum(!is.na(background))
#                   }
#                }
#                
#                # Subset dn:
#                if(!is.na(nGenesDn[iGeneSet,iContrast])) {
#                   if(nGenesDn[iGeneSet,iContrast] != 0) {
#                      value <- gsStatsDn[iGeneSet,iContrast]  
#                      background <- gsStatsDnPerm[[iContrast]][iGeneSet,]
#                      pValuesDn[iGeneSet] <- sum(background >= value, na.rm=TRUE)/sum(!is.na(background))
#                   }
#                }
#             }
         }
         
         
         #************************  
         # Wilcoxon:
         #************************
         if(statMethod == "wilcoxon") {
            
            # "Directional":
            if(statType == "t") {
               value <- gsStatsAll[iGeneSet,iContrast]  
               background <- gsStatsAllPerm[[iContrast]][iGeneSet,]
               pValuesAllUp[iGeneSet] <- sum(background >= value, na.rm=TRUE)/sum(!is.na(background)) # all up
               pValuesAllDn[iGeneSet] <- sum(background <= value, na.rm=TRUE)/sum(!is.na(background)) # all dn
            } else if(statType == "p-signed") {
               valueTestUp <- gsStatsAllTestUp[iGeneSet,iContrast]
               valueTestDn <- gsStatsAllTestDn[iGeneSet,iContrast]
               backgroundTestUp <- gsStatsAllTestUpPerm[[iContrast]][iGeneSet,]
               backgroundTestDn <- gsStatsAllTestDnPerm[[iContrast]][iGeneSet,]
               pValuesAllUp[iGeneSet] <- sum(backgroundTestUp <= valueTestUp, na.rm=TRUE)/sum(!is.na(backgroundTestUp)) # all up
               pValuesAllDn[iGeneSet] <- sum(backgroundTestDn <= valueTestDn, na.rm=TRUE)/sum(!is.na(backgroundTestDn)) # all dn
            }
            
            # "Mixed":
            value <- gsStatsAbs[iGeneSet,iContrast]  
            background <- gsStatsAbsPerm[[iContrast]][iGeneSet,]
            if(statType %in% c("p","p-signed")) {
               pValuesAbs[iGeneSet] <- sum(background <= value, na.rm=TRUE)/sum(!is.na(background))
            } else {
               pValuesAbs[iGeneSet] <- sum(background >= value, na.rm=TRUE)/sum(!is.na(background))
            }
            
#             # Subset up:
#             if(!is.na(nGenesUp[iGeneSet,iContrast])) {
#                if(nGenesUp[iGeneSet,iContrast] != 0) {
#                   value <- gsStatsUp[iGeneSet,iContrast]  
#                   background <- gsStatsUpPerm[[iContrast]][iGeneSet,]
#                   if(statType == "p-signed") {
#                      pValuesUp[iGeneSet] <- sum(background <= value, na.rm=TRUE)/sum(!is.na(background))
#                   } else if(statType %in% c("t","F-signed")) {
#                      pValuesUp[iGeneSet] <- sum(background >= value, na.rm=TRUE)/sum(!is.na(background))
#                   }
#                }
#             }
#             
#             # Subset dn:
#             if(!is.na(nGenesDn[iGeneSet,iContrast])) {
#                if(nGenesDn[iGeneSet,iContrast] != 0) {
#                   value <- gsStatsDn[iGeneSet,iContrast]  
#                   background <- gsStatsDnPerm[[iContrast]][iGeneSet,]
#                   if(statType == "p-signed") {
#                      pValuesDn[iGeneSet] <- sum(background <= value, na.rm=TRUE)/sum(!is.na(background))
#                   } else if(statType %in% c("t","F-signed")) {
#                      pValuesDn[iGeneSet] <- sum(background >= value, na.rm=TRUE)/sum(!is.na(background))
#                   }
#                }  
#             }
         }
         
         #************************
         # Mean:
         # Median:
         # Sum:
         #************************
         if(statMethod %in% c("mean","median","sum")) {
            
            # "Directional":
            if(statType == "t") {
               value <- gsStatsAll[iGeneSet,iContrast]
               background <- gsStatsAllPerm[[iContrast]][iGeneSet,]
               pValuesAllUp[iGeneSet] <- sum(background >= value, na.rm=TRUE)/sum(!is.na(background)) # all up
               pValuesAllDn[iGeneSet] <- sum(background <= value, na.rm=TRUE)/sum(!is.na(background)) # all dn
            } else if(statType == "p-signed") {
               valueTestUp <- gsStatsAllTestUp[iGeneSet,iContrast]
               valueTestDn <- gsStatsAllTestDn[iGeneSet,iContrast]
               backgroundTestUp <- gsStatsAllTestUpPerm[[iContrast]][iGeneSet,]
               backgroundTestDn <- gsStatsAllTestDnPerm[[iContrast]][iGeneSet,]
               pValuesAllUp[iGeneSet] <- sum(backgroundTestUp <= valueTestUp, na.rm=TRUE)/sum(!is.na(backgroundTestUp)) # all up
               pValuesAllDn[iGeneSet] <- sum(backgroundTestDn <= valueTestDn, na.rm=TRUE)/sum(!is.na(backgroundTestDn)) # all dn
            }
            
            # "Mixed":
            value <- gsStatsAbs[iGeneSet,iContrast]  
            background <- gsStatsAbsPerm[[iContrast]][iGeneSet,]
            if(statType %in% c("p","p-signed")) {
               pValuesAbs[iGeneSet] <- sum(background <= value, na.rm=TRUE)/sum(!is.na(background))
            } else {
               pValuesAbs[iGeneSet] <- sum(background >= value, na.rm=TRUE)/sum(!is.na(background))
            }
            
#             # Subset up:
#             if(!is.na(nGenesUp[iGeneSet,iContrast])) {
#                if(nGenesUp[iGeneSet,iContrast] != 0) {
#                   value <- gsStatsUp[iGeneSet,iContrast]  
#                   background <- gsStatsUpPerm[[iContrast]][iGeneSet,]
#                   if(statType == "p-signed") {
#                      pValuesUp[iGeneSet] <- sum(background <= value, na.rm=TRUE)/sum(!is.na(background))
#                   } else if(statType %in% c("t","F-signed")){
#                      pValuesUp[iGeneSet] <- sum(background >= value, na.rm=TRUE)/sum(!is.na(background))
#                   }
#                }
#             }
#             
#             # Subset dn:
#             if(!is.na(nGenesDn[iGeneSet,iContrast])) {
#                if(nGenesDn[iGeneSet,iContrast] != 0) {
#                   value <- gsStatsDn[iGeneSet,iContrast]  
#                   background <- gsStatsDnPerm[[iContrast]][iGeneSet,]
#                   if(statType == "p-signed") {
#                      pValuesDn[iGeneSet] <- sum(background <= value, na.rm=TRUE)/sum(!is.na(background))
#                   } else if(statType %in% c("t","F-signed")){
#                      pValuesDn[iGeneSet] <- sum(background >= value, na.rm=TRUE)/sum(!is.na(background))
#                   }
#                }
#             }
         }
         
         #************************   
         # Maxmean:
         #************************
         if(statMethod == "maxmean") {
            
            # "Mixed":
            value <- gsStatsAbs[iGeneSet,iContrast]  
            background <- gsStatsAbsPerm[[iContrast]][iGeneSet,]
            pValuesAbs[iGeneSet] <- sum(background >= value, na.rm=TRUE)/sum(!is.na(background)) 
         }
         
         #************************   
         # GSEA:
         #************************
         if(statMethod == "gsea") {
            value <- gsStatsAll[iGeneSet,iContrast]   
            background <- gsStatsAllPerm[[iContrast]][iGeneSet,]
            if(value < 0) {
               backgroundNeg <- background[background < 0]
               pValuesAllDn[iGeneSet] <- sum(backgroundNeg <= value, na.rm=TRUE)/sum(!is.na(backgroundNeg)) # all
               pValuesAllUp[iGeneSet] <- NA
            } else {
               backgroundPos <- background[background > 0]
               pValuesAllUp[iGeneSet] <- sum(backgroundPos >= value, na.rm=TRUE)/sum(!is.na(backgroundPos)) # all
               pValuesAllDn[iGeneSet] <- NA
            }
            
         }
         
         #************************  
         # Page:
         #************************
         if(statMethod == "page") {
            value <- gsStatsAll[iGeneSet,iContrast]    
            background <- gsStatsAllPerm[[iContrast]][iGeneSet,]
            pValuesAllUp[iGeneSet] <- sum(background >= value, na.rm=TRUE)/sum(!is.na(background)) # all up
            pValuesAllDn[iGeneSet] <- sum(background <= value, na.rm=TRUE)/sum(!is.na(background)) # all dn
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