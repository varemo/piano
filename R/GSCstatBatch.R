GSCstatBatch <- function(statistics, statType, gsc, method, signMethod, gseaParam, signs) {
   
   # Create result list object:
   res <- list()
   
   res$statsAll <- vector()
   res$statsAllTestUp <- vector()
   res$statsAllTestDn <- vector()
   res$statsAbs <- vector()
   res$statsUp <- vector()
   res$statsDn <- vector()
   
   res$nGenes <- vector()
   res$nGenesUp <- vector()
   res$nGenesDn <- vector()
   
   res$pValuesAll <- vector()
   res$pValuesAllUp <- vector()
   res$pValuesAllDn <- vector()
   res$pValuesAbs <- vector()
   res$pValuesUp <- vector()
   res$pValuesDn <- vector()
   
   # Run separately for each contrast/comparison (statistics column):
   for(iContrast in 1:ncol(statistics)) {
      
      # Preallocate:
      gsStatsAll <- rep(NA,length(gsc))
      gsStatsAllTestUp <- rep(NA,length(gsc))
      gsStatsAllTestDn <- rep(NA,length(gsc))
      gsStatsAbs <- rep(NA,length(gsc))
      gsStatsUp <- rep(NA,length(gsc))
      gsStatsDn <- rep(NA,length(gsc))
      nGenes <- rep(0,length(gsc))
      nGenesUp <- rep(0,length(gsc))
      nGenesDn <- rep(0,length(gsc))
      
      pValuesAll <- rep(NA,length(gsc))
      pValuesAllUp <- rep(NA,length(gsc))
      pValuesAllDn <- rep(NA,length(gsc))
      pValuesAbs <- rep(NA,length(gsc))
      pValuesUp <- rep(NA,length(gsc))
      pValuesDn <- rep(NA,length(gsc))
      
      # Get statistics for current contrast:
      statsContrast <- statistics[,iContrast]
      if(statType %in% c("p-signed","F-signed")) {
         signsContrast <- signs[,iContrast]
      }
      if(statType == "t") {
         signsContrast <- sign(statsContrast)
      }
      
      # If possible, calculate "directional p-values":
      if(method %in% c("stouffer","reporter","tailStrength","wilcoxon","wilcoxon_fast","mean","median","sum") & statType == "p-signed") {
         statsContrastTestUp <- abs(statsContrast)
         statsContrastTestUp[signsContrast > 0] <- statsContrastTestUp[signsContrast > 0]/2
         statsContrastTestUp[signsContrast < 0] <- 1 - statsContrastTestUp[signsContrast < 0]/2
         statsContrastTestUp[signsContrast == 0] <- NA
         statsContrastTestDn <- 1 - statsContrastTestUp
         
         # Fix zeros and ones in directional p-values:
         statsContrastTestUp[statsContrastTestUp < 1e-100] <- 1e-100
         statsContrastTestUp[statsContrastTestUp > 0.9999999999999999] <- 0.9999999999999999
         statsContrastTestDn[statsContrastTestDn < 1e-100] <- 1e-100
         statsContrastTestDn[statsContrastTestDn > 0.9999999999999999] <- 0.9999999999999999
      }
      
      # Loop through the gene-sets and calculate gene-set statistics:
      for(iGeneSet in 1:length(gsc)) {
         
         #************************
         # Fisher:
         #************************
         if(method == "fisher") {
            
            # Get statistics for gene-set:
            indGenesInSet <- names(statsContrast) %in% gsc[[iGeneSet]]
            statsGenesInSet <- statsContrast[indGenesInSet]
            nGenes[iGeneSet] <- length(statsGenesInSet)
            if(nGenes[iGeneSet] > 0) {
               
               # Genes up/dn (partly redundant code added later as a fix):
               if(exists("signsContrast")) {
                  nGenesUp[iGeneSet] <- sum(signsContrast[indGenesInSet]>0)
                  nGenesDn[iGeneSet] <- sum(signsContrast[indGenesInSet]<0)
               }
               
               # "Mixed" gene-set statistic:
               gsStatsAbs[iGeneSet] <- calcGeneSetStat(abs(statsGenesInSet), method)
               
               # "Subset-up" gene-set statistic:
               if(statType == "p-signed" & signMethod != "sampleperm") {
                  signsGenesInSet <- signsContrast[indGenesInSet]
                  sel <- statsGenesInSet[signsGenesInSet > 0]
                  nGenesUp[iGeneSet] <- length(sel)
                  if(nGenesUp[iGeneSet] > 0) {
                     gsStatsUp[iGeneSet] <- calcGeneSetStat(sel, method)
                  }
                  
                  # "Subset-dn" gene-set statistic:
                  sel <- abs(statsGenesInSet[signsGenesInSet < 0])
                  nGenesDn[iGeneSet] <- length(sel)
                  if(nGenesDn[iGeneSet] > 0) {
                     gsStatsDn[iGeneSet] <- calcGeneSetStat(sel, method)  
                  }
               } 
            }
            
            
         #************************
         # Stouffer:
         # Reporter:
         # Tail strength:
         #************************
         } else if(method %in% c("stouffer","reporter","tailStrength")) {
            
            # Get statistics for gene-set:
            indGenesInSet <- names(statsContrast) %in% gsc[[iGeneSet]]
            statsGenesInSet <- statsContrast[indGenesInSet]
            nGenes[iGeneSet] <- length(statsGenesInSet)
            if(nGenes[iGeneSet] > 0) {
               
               # Genes up/dn (partly redundant code added later as a fix):
               if(exists("signsContrast")) {
                  nGenesUp[iGeneSet] <- sum(signsContrast[indGenesInSet]>0)
                  nGenesDn[iGeneSet] <- sum(signsContrast[indGenesInSet]<0)
               }
               
               # "Directional" gene-set statistics:
               if(statType == "p-signed") {
                  statsGenesInSetTestUp <- statsContrastTestUp[indGenesInSet]
                  statsGenesInSetTestDn <- statsContrastTestDn[indGenesInSet]
                  statsGenesInSetTestUp <- statsGenesInSetTestUp[!is.na(statsGenesInSetTestUp)]
                  statsGenesInSetTestDn <- statsGenesInSetTestDn[!is.na(statsGenesInSetTestDn)]
                  gsStatsAllTestUp[iGeneSet] <- calcGeneSetStat(statsGenesInSetTestUp, method)
                  gsStatsAllTestDn[iGeneSet] <- calcGeneSetStat(statsGenesInSetTestDn, method)
               }
               
               # "Mixed" gene-set statistic:
               gsStatsAbs[iGeneSet] <- calcGeneSetStat(abs(statsGenesInSet), method)
               
               # "Subset-up" gene-set statistic:
               if(statType == "p-signed" & signMethod != "sampleperm") {
                  signsGenesInSet <- signsContrast[indGenesInSet]
                  sel <- statsGenesInSet[signsGenesInSet > 0]
                  nGenesUp[iGeneSet] <- length(sel)
                  if(nGenesUp[iGeneSet] > 0) {
                     gsStatsUp[iGeneSet] <- calcGeneSetStat(sel, method)
                  }
                  
                  # "Subset-dn" gene-set statistic:
                  sel <- abs(statsGenesInSet[signsGenesInSet < 0])
                  nGenesDn[iGeneSet] <- length(sel)
                  if(nGenesDn[iGeneSet] > 0) {
                     gsStatsDn[iGeneSet] <- calcGeneSetStat(sel, method)  
                  }
               }
            }
         
         
         #************************
         # Wilcoxon:
         #************************
         } else if(method == "wilcoxon") {
            
            # Get statistics for gene-set:
            indGenesInSet <- names(statsContrast) %in% gsc[[iGeneSet]]
            statsGenesInSet <- statsContrast[indGenesInSet]
            statsNotInSet <- statsContrast[!indGenesInSet]
            if(statType %in% c("p-signed","F-signed")) {
               signsGenesInSet <- signsContrast[indGenesInSet]
               signsNotInSet <- signsContrast[!indGenesInSet]
            }
            nGenes[iGeneSet] <- length(statsGenesInSet)
            if(nGenes[iGeneSet] > 0) {
              
               # Genes up/dn (partly redundant code added later as a fix):
               if(exists("signsContrast")) {
                  nGenesUp[iGeneSet] <- sum(signsContrast[indGenesInSet]>0)
                  nGenesDn[iGeneSet] <- sum(signsContrast[indGenesInSet]<0)
               }
               
               # "Directional" gene-set statistics:
               if(statType == "p-signed") {
                  statsGenesInSetTestUp <- statsContrastTestUp[indGenesInSet]
                  statsGenesInSetTestDn <- statsContrastTestDn[indGenesInSet]
                  statsGenesInSetTestUp <- statsGenesInSetTestUp[!is.na(statsGenesInSetTestUp)]
                  statsGenesInSetTestDn <- statsGenesInSetTestDn[!is.na(statsGenesInSetTestDn)]
                  statsNotInSetTestUp <- statsContrastTestUp[!indGenesInSet]
                  statsNotInSetTestDn <- statsContrastTestDn[!indGenesInSet]
                  statsNotInSetTestUp <- statsNotInSetTestUp[!is.na(statsNotInSetTestUp)]
                  statsNotInSetTestDn <- statsNotInSetTestDn[!is.na(statsNotInSetTestDn)]
                  tmp <- calcGeneSetStat(statsGenesInSetTestUp, "wilcoxon_less", statsNotInSetTestUp)
                  gsStatsAllTestUp[iGeneSet] <- tmp[1]
                  pValuesAllUp[iGeneSet] <- tmp[2]
                  tmp <- calcGeneSetStat(statsGenesInSetTestDn, "wilcoxon_less", statsNotInSetTestDn)
                  gsStatsAllTestDn[iGeneSet] <- tmp[1]
                  pValuesAllDn[iGeneSet] <- tmp[2]
               } else if(statType == "t") {
                  tmp <- calcGeneSetStat(statsGenesInSet, "wilcoxon_two.sided", statsNotInSet)
                  gsStatsAll[iGeneSet] <- tmp[1]
                  #pValuesAll[iGeneSet] <- tmp[2]
                  tmp <- calcGeneSetStat(statsGenesInSet, "wilcoxon_greater", statsNotInSet)
                  pValuesAllUp[iGeneSet] <- tmp[2]
                  tmp <- calcGeneSetStat(statsGenesInSet, "wilcoxon_less", statsNotInSet)
                  pValuesAllDn[iGeneSet] <- tmp[2]   
               }
               
               
               # "Mixed" gene-set statistic:
               if(statType %in% c("p","p-signed")) {
                  tmp <- calcGeneSetStat(abs(statsGenesInSet), "wilcoxon_less", abs(statsNotInSet))
                  gsStatsAbs[iGeneSet] <- tmp[1]
                  pValuesAbs[iGeneSet] <- tmp[2]
               } else if(statType %in% c("t","F","F-signed")) {
                  tmp <- calcGeneSetStat(abs(statsGenesInSet), "wilcoxon_greater", abs(statsNotInSet))
                  gsStatsAbs[iGeneSet] <- tmp[1]
                  pValuesAbs[iGeneSet] <- tmp[2]
               }
               
               # "Subset-up" gene-set statistic:
               if(statType %in% c("t","p-signed","F-signed") & signMethod != "sampleperm") {
                  if(statType == "t") {
                     sel <- statsGenesInSet[statsGenesInSet > 0]
                  } else {
                     sel <- statsGenesInSet[signsGenesInSet > 0]
                  }
                  nGenesUp[iGeneSet] <- length(sel)
                  if(nGenesUp[iGeneSet] > 0) {
                     if(statType == "p-signed") {
                        tmp <- calcGeneSetStat(sel, "wilcoxon_less", statsNotInSet[signsNotInSet > 0])
                        gsStatsUp[iGeneSet] <- tmp[1]
                        pValuesUp[iGeneSet] <- tmp[2]
                     } else if(statType == "t") {
                        tmp <- calcGeneSetStat(sel, "wilcoxon_greater",  statsNotInSet[statsNotInSet > 0])
                        gsStatsUp[iGeneSet] <- tmp[1]
                        pValuesUp[iGeneSet] <- tmp[2] 
                     } else if(statType == "F-signed") {
                        tmp <- calcGeneSetStat(sel, "wilcoxon_greater",  statsNotInSet[signsNotInSet > 0])
                        gsStatsUp[iGeneSet] <- tmp[1]
                        pValuesUp[iGeneSet] <- tmp[2] 
                     }
                  }
                  
                  # "Subset-dn" gene-set statistic:
                  if(statType == "t") {
                     sel <- abs(statsGenesInSet[statsGenesInSet < 0])
                  } else {
                     sel <- abs(statsGenesInSet[signsGenesInSet < 0])
                  }
                  nGenesDn[iGeneSet] <- length(sel)
                  if(nGenesDn[iGeneSet] > 0) {
                     if(statType == "p-signed") {
                        tmp <- calcGeneSetStat(sel, "wilcoxon_less", abs(statsNotInSet[signsNotInSet < 0]))
                        gsStatsDn[iGeneSet] <- tmp[1]
                        pValuesDn[iGeneSet] <- tmp[2]
                     } else if(statType == "t") {
                        tmp <- calcGeneSetStat(sel, "wilcoxon_greater",  abs(statsNotInSet[statsNotInSet < 0]))
                        gsStatsDn[iGeneSet] <- tmp[1]
                        pValuesDn[iGeneSet] <- tmp[2] 
                     } else if(statType == "F-signed") {
                        tmp <- calcGeneSetStat(sel, "wilcoxon_greater",  abs(statsNotInSet[signsNotInSet < 0]))
                        gsStatsDn[iGeneSet] <- tmp[1]
                        pValuesDn[iGeneSet] <- tmp[2] 
                     }
                  }
               }
            }
            
            
         #************************
         # Wilcoxon, fast (only used by GSCstatSamplePerm):
         #************************
         } else if(method == "wilcoxon_fast") {
            
            # Get statistics for gene-set:
            indGenesInSet <- names(statsContrast) %in% gsc[[iGeneSet]]
            statsGenesInSet <- statsContrast[indGenesInSet]
            statsNotInSet <- statsContrast[!indGenesInSet]
            if(statType %in% c("p-signed","F-signed")) {
               signsGenesInSet <- signsContrast[indGenesInSet]
               signsNotInSet <- signsContrast[!indGenesInSet]
            }
            nGenes[iGeneSet] <- length(statsGenesInSet)
            if(nGenes[iGeneSet] > 0) {
               
               # Genes up/dn (partly redundant code added later as a fix):
               if(exists("signsContrast")) {
                  nGenesUp[iGeneSet] <- sum(signsContrast[indGenesInSet]>0)
                  nGenesDn[iGeneSet] <- sum(signsContrast[indGenesInSet]<0)
               }
               
               # "Directional" gene-set statistics:
               if(statType == "p-signed") {
                  statsGenesInSetTestUp <- statsContrastTestUp[indGenesInSet]
                  statsGenesInSetTestDn <- statsContrastTestDn[indGenesInSet]
                  statsGenesInSetTestUp <- statsGenesInSetTestUp[!is.na(statsGenesInSetTestUp)]
                  statsGenesInSetTestDn <- statsGenesInSetTestDn[!is.na(statsGenesInSetTestDn)]
                  statsNotInSetTestUp <- statsContrastTestUp[!indGenesInSet]
                  statsNotInSetTestDn <- statsContrastTestDn[!indGenesInSet]
                  statsNotInSetTestUp <- statsNotInSetTestUp[!is.na(statsNotInSetTestUp)]
                  statsNotInSetTestDn <- statsNotInSetTestDn[!is.na(statsNotInSetTestDn)]
                  gsStatsAllTestUp[iGeneSet] <- calcGeneSetStat(statsGenesInSetTestUp, "wilcoxon_fast", statsNotInSetTestUp)                              
                  gsStatsAllTestDn[iGeneSet] <- calcGeneSetStat(statsGenesInSetTestDn, "wilcoxon_fast", statsNotInSetTestDn)                                    
               } else if(statType == "t") {
                  gsStatsAll[iGeneSet] <- calcGeneSetStat(statsGenesInSet, "wilcoxon_fast", statsNotInSet)                   
               }
               
               
               # "Mixed" gene-set statistic:
               if(statType %in% c("p","p-signed")) {
                  gsStatsAbs[iGeneSet] <- calcGeneSetStat(abs(statsGenesInSet), "wilcoxon_fast", abs(statsNotInSet))                                     
               } else if(statType %in% c("t","F","F-signed")) {
                  gsStatsAbs[iGeneSet] <- calcGeneSetStat(abs(statsGenesInSet), "wilcoxon_fast", abs(statsNotInSet))                              
               }
               
               # "Subset-up" gene-set statistic:
               if(statType %in% c("t","p-signed","F-signed") & signMethod != "sampleperm") {
                  if(statType == "t") {
                     sel <- statsGenesInSet[statsGenesInSet > 0]
                  } else {
                     sel <- statsGenesInSet[signsGenesInSet > 0]
                  }
                  nGenesUp[iGeneSet] <- length(sel)
                  if(nGenesUp[iGeneSet] > 0) {
                     if(statType == "p-signed") {
                        gsStatsUp[iGeneSet] <- calcGeneSetStat(sel, "wilcoxon_fast", statsNotInSet[signsNotInSet > 0])                     
                     } else if(statType == "t") {
                        gsStatsUp[iGeneSet] <- calcGeneSetStat(sel, "wilcoxon_fast",  statsNotInSet[statsNotInSet > 0])                                                
                     } else if(statType == "F-signed") {
                        gsStatsUp[iGeneSet] <- calcGeneSetStat(sel, "wilcoxon_fast",  statsNotInSet[signsNotInSet > 0])                        
                     }
                  }
                  
                  # "Subset-dn" gene-set statistic:
                  if(statType == "t") {
                     sel <- abs(statsGenesInSet[statsGenesInSet < 0])
                  } else {
                     sel <- abs(statsGenesInSet[signsGenesInSet < 0])
                  }
                  nGenesDn[iGeneSet] <- length(sel)
                  if(nGenesDn[iGeneSet] > 0) {
                     if(statType == "p-signed") {
                        gsStatsDn[iGeneSet] <- calcGeneSetStat(sel, "wilcoxon_fast", abs(statsNotInSet[signsNotInSet < 0]))                                              
                     } else if(statType == "t") {
                        gsStatsDn[iGeneSet] <- calcGeneSetStat(sel, "wilcoxon_fast",  abs(statsNotInSet[statsNotInSet < 0]))                                                 
                     } else if(statType == "F-signed") {
                        gsStatsDn[iGeneSet] <- calcGeneSetStat(sel, "wilcoxon_fast",  abs(statsNotInSet[signsNotInSet < 0]))                        
                     }
                  }
               }
            }
            
         
         #************************
         # Mean:
         # Median:
         # Sum:
         #************************
         } else if(method %in% c("mean","median","sum")) {

            # Get statistics for gene-set:
            indGenesInSet <- names(statsContrast) %in% gsc[[iGeneSet]]
            statsGenesInSet <- statsContrast[indGenesInSet]
            if(statType %in% c("p-signed","F-signed")) {
               signsGenesInSet <- signsContrast[indGenesInSet]
            }
            nGenes[iGeneSet] <- length(statsGenesInSet)
            if(nGenes[iGeneSet] > 0) {
               
               # Genes up/dn (partly redundant code added later as a fix):
               if(exists("signsContrast")) {
                  nGenesUp[iGeneSet] <- sum(signsContrast[indGenesInSet]>0)
                  nGenesDn[iGeneSet] <- sum(signsContrast[indGenesInSet]<0)
               }
               
               # "Directional" gene-set statistics:
               if(statType == "p-signed") {
                  statsGenesInSetTestUp <- statsContrastTestUp[indGenesInSet]
                  statsGenesInSetTestDn <- statsContrastTestDn[indGenesInSet]
                  statsGenesInSetTestUp <- statsGenesInSetTestUp[!is.na(statsGenesInSetTestUp)]
                  statsGenesInSetTestDn <- statsGenesInSetTestDn[!is.na(statsGenesInSetTestDn)]
                  gsStatsAllTestUp[iGeneSet] <- calcGeneSetStat(statsGenesInSetTestUp, method)
                  gsStatsAllTestDn[iGeneSet] <- calcGeneSetStat(statsGenesInSetTestDn, method)
               } else if(statType == "t") {
                  gsStatsAll[iGeneSet] <- calcGeneSetStat(statsGenesInSet, method)  
               }
               
               # "Mixed" gene-set statistic:
               gsStatsAbs[iGeneSet] <- calcGeneSetStat(abs(statsGenesInSet), method)
               
               # "Subset-up" gene-set statistic:
               if(statType %in% c("p-signed","t","F-signed") & signMethod != "sampleperm") {
                  if(statType == "t") {
                     sel <- statsGenesInSet[statsGenesInSet > 0]
                  } else {
                     sel <- statsGenesInSet[signsGenesInSet > 0]
                  }
                  nGenesUp[iGeneSet] <- length(sel)
                  if(nGenesUp[iGeneSet] > 0) {
                     gsStatsUp[iGeneSet] <- calcGeneSetStat(sel, method)
                  }
                  
                  # "Subset-dn" gene-set statistic:
                  if(statType == "t") {
                     sel <- abs(statsGenesInSet[statsGenesInSet < 0])
                  } else {
                     sel <- abs(statsGenesInSet[signsGenesInSet < 0])
                  }
                  nGenesDn[iGeneSet] <- length(sel)
                  if(nGenesDn[iGeneSet] > 0) {
                     gsStatsDn[iGeneSet] <- calcGeneSetStat(sel, method)  
                  }
               }
            }
            
         #************************   
         # Maxmean:
         #************************
         } else if(method == "maxmean") {
            
            # Get statistics for gene-set:
            indGenesInSet <- names(statsContrast) %in% gsc[[iGeneSet]]
            statsGenesInSet <- statsContrast[indGenesInSet]
            nGenes[iGeneSet] <- length(statsGenesInSet)
            if(nGenes[iGeneSet] > 0) {
               
               # Genes up/dn (partly redundant code added later as a fix):
               if(exists("signsContrast")) {
                  nGenesUp[iGeneSet] <- sum(signsContrast[indGenesInSet]>0)
                  nGenesDn[iGeneSet] <- sum(signsContrast[indGenesInSet]<0)
               }
               
               # "Mixed" gene-set statistic:
               gsStatsAbs[iGeneSet] <- calcGeneSetStat(statsGenesInSet, method)
            }
            
         #************************   
         # GSEA:
         #************************
         } else if(method == "gsea") {
            
            # Get statistics for gene-set:
            statsContrastSorted <- sort(statsContrast,decreasing=TRUE)
            indGenesInSet <- which(names(statsContrastSorted) %in% gsc[[iGeneSet]])
            nGenes[iGeneSet] <- length(indGenesInSet)
            if(nGenes[iGeneSet] > 0) {
               
               # Genes up/dn (partly redundant code added later as a fix):
               if(exists("signsContrast")) {
                  #nGenesUp[iGeneSet] <- sum(signsContrast[indGenesInSet]>0)
                  #nGenesDn[iGeneSet] <- sum(signsContrast[indGenesInSet]<0)
                  nGenesUp[iGeneSet] <- sum(statsContrastSorted[indGenesInSet]>0)
                  nGenesDn[iGeneSet] <- sum(statsContrastSorted[indGenesInSet]<0)
               }
               
               # "Directional" gene-set statistic:
               gsStatsAll[iGeneSet] <- calcGeneSetStat(indGenesInSet, method, statsContrastSorted, gseaParam)            
            }
            
            
         #************************  
         # Page:
         #************************
         } else if(method == "page") {
            
            # Get statistic for gene-set:
            indGenesInSet <- names(statsContrast) %in% gsc[[iGeneSet]]
            statsGenesInSet <- statsContrast[indGenesInSet]
            nGenes[iGeneSet] <- length(statsGenesInSet)
            if(nGenes[iGeneSet] > 0) {
               
               # Genes up/dn (partly redundant code added later as a fix):
               if(exists("signsContrast")) {
                  nGenesUp[iGeneSet] <- sum(signsContrast[indGenesInSet]>0)
                  nGenesDn[iGeneSet] <- sum(signsContrast[indGenesInSet]<0)
               }
               
               # "Directional" gene-set statistics:
               gsStatsAll[iGeneSet] <- calcGeneSetStat(statsGenesInSet, method, statsContrast)
            }
         }
      }
      
      # Save results:
      res$statsAll <- cbind(res$statsAll, gsStatsAll)
      res$statsAllTestUp <- cbind(res$statsAllTestUp, gsStatsAllTestUp)
      res$statsAllTestDn <- cbind(res$statsAllTestDn, gsStatsAllTestDn)
      res$statsAbs <- cbind(res$statsAbs, gsStatsAbs)
      res$statsUp <- cbind(res$statsUp, gsStatsUp)
      res$statsDn <- cbind(res$statsDn, gsStatsDn)
      
      res$nGenes <- cbind(res$nGenes, nGenes)
      res$nGenesUp <- cbind(res$nGenesUp, nGenesUp)
      res$nGenesDn <- cbind(res$nGenesDn, nGenesDn)
      
      res$pValuesAll <- cbind(res$pValuesAll, pValuesAll)
      res$pValuesAllUp <- cbind(res$pValuesAllUp, pValuesAllUp)
      res$pValuesAllDn <- cbind(res$pValuesAllDn, pValuesAllDn)
      res$pValuesAbs <- cbind(res$pValuesAbs, pValuesAbs)
      res$pValuesUp <- cbind(res$pValuesUp, pValuesUp)
      res$pValuesDn <- cbind(res$pValuesDn, pValuesDn)
   }
   
   if(method == "fisher") res$statName = "X2"
   else if(method == "stouffer") res$statName = "Z"
   else if(method == "tailStrength") res$statName = "TS"
   else if(method == "reporter") res$statName = "Z" # This is the raw Z-score.
   else if(method == "wilcoxon") res$statName = "W"
   else if(method == "wilcoxon_fast") res$statName = "W"
   else if(method == "mean") res$statName = "mean"
   else if(method == "median") res$statName = "median"
   else if(method == "sum") res$statName = "sum"
   else if(method == "maxmean") res$statName = "maxmean"
   else if(method == "gsea") res$statName = "ES"
   else if(method == "page") res$statName = "Z"
   return(res)
   
}