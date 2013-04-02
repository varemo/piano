GSCstatGenePerm <- function(statistics, signs, gsc, statType, method, nGenes, nGenesUp, nGenesDn, nPerm, gseaParam) {

   # Create result list object:
   res <- list()
   res$gsStatsAllPerm <- list()
   res$gsStatsAllTestUpPerm <- list()
   res$gsStatsAllTestDnPerm <- list()
   res$gsStatsAbsPerm <- list()
   res$gsStatsUpPerm <- list()
   res$gsStatsDnPerm <- list()
   
   # Run separately for each contrast/comparison (statistics column):
   for(iContrast in 1:ncol(statistics)) {
      
      # Get statistics for current contrast:
      statsContrast <- statistics[,iContrast] # All
      if(statType %in% c("p-signed","F-signed")) {
         signsContrast <- signs[,iContrast]
         statsContrastUp <- statsContrast[signsContrast > 0] # Subset up
         statsContrastDn <- abs(statsContrast[signsContrast < 0]) # Subset dn
      } else {
         statsContrastUp <- statsContrast[statsContrast > 0] # Subset up
         statsContrastDn <- abs(statsContrast[statsContrast < 0]) # Subset dn
      }
      if(method == "gsea") {
         statsContrastSorted <- sort(statistics[,iContrast],decreasing=TRUE)
      }
      
      # If possible, calculate "directional p-values":
      if(method %in% c("stouffer","reporter","tailStrength","wilcoxon","mean","median","sum") & statType == "p-signed") {
         statsContrastTestUp <- abs(statsContrast)
         statsContrastTestUp[signsContrast > 0] <- statsContrastTestUp[signsContrast > 0]/2
         statsContrastTestUp[signsContrast < 0] <- 1 - statsContrastTestUp[signsContrast < 0]/2
         statsContrastTestUp[signsContrast == 0] <- NA
         statsContrastTestUp <- statsContrastTestUp[!is.na(statsContrastTestUp)]
         statsContrastTestDn <- 1 - statsContrastTestUp
         
         # Fix zeros and ones in directional p-values:
         statsContrastTestUp[statsContrastTestUp < 1e-100] <- 1e-100
         statsContrastTestUp[statsContrastTestUp > 0.9999999999999999] <- 0.9999999999999999
         statsContrastTestDn[statsContrastTestDn < 1e-100] <- 1e-100
         statsContrastTestDn[statsContrastTestDn > 0.9999999999999999] <- 0.9999999999999999
      }
      
      # Get gene set sizes for current contrast:
      gsSizes <- unique(nGenes[,iContrast])
      gsSizesUp <- unique(nGenesUp[,iContrast])
      gsSizesDn <- unique(nGenesDn[,iContrast])
      
      # Preallocate:
      gsStatsAllMatrix <- matrix(nrow=length(gsSizes), ncol=nPerm)
      gsStatsAllTestUpMatrix <- matrix(nrow=length(gsSizes), ncol=nPerm)
      gsStatsAllTestDnMatrix <- matrix(nrow=length(gsSizes), ncol=nPerm)
      gsStatsAbsMatrix <- matrix(nrow=length(gsSizes), ncol=nPerm)
      gsStatsUpMatrix <- matrix(nrow=length(gsSizesUp), ncol=nPerm)
      gsStatsDnMatrix <- matrix(nrow=length(gsSizesDn), ncol=nPerm)
      
      
      # For each full gene set size (Directional and mixed):
      for(iGeneSet in 1:length(gsSizes)) {
         gsStatsAll <- rep(NA,nPerm)
         gsStatsAllTestUp <- rep(NA,nPerm)
         gsStatsAllTestDn <- rep(NA,nPerm)
         gsStatsAbs <- rep(NA,nPerm)
         nGenesInSet <- gsSizes[iGeneSet]
         
         # Loop nPerm times:
         for(iPerm in 1:nPerm) {
            
            #************************  
            # Fisher:
            #************************
            if(method == "fisher") {
            
               # "Mixed" gene-set statistics:
               #rInd <- ceiling(runif(nGenesInSet,0,length(statsContrast)))
               rInd <- sample(1:length(statsContrast),nGenesInSet)
               randStats <- statsContrast[rInd]
               gsStatsAbs[iPerm] <- calcGeneSetStat(abs(randStats), method)
            }
              
            
            #************************
            # Stouffer:
            # Reporter:
            # Tail strength:
            #************************
            if(method %in% c("stouffer","reporter","tailStrength")) {
               
               # Select random set:
               #rInd <- ceiling(runif(nGenesInSet,0,length(statsContrast)))
               rInd <- sample(1:length(statsContrast),nGenesInSet)
               randStats <- statsContrast[rInd]
               
               # "Directional" gene-set statistics:
               if(statType == "p-signed") {
                  rInd <- sample(1:length(statsContrastTestUp),nGenesInSet)
                  randStatsTestUp <- statsContrastTestUp[rInd]
                  randStatsTestDn <- statsContrastTestDn[rInd]                  
                  gsStatsAllTestUp[iPerm] <- calcGeneSetStat(randStatsTestUp, method)                  
                  gsStatsAllTestDn[iPerm] <- calcGeneSetStat(randStatsTestDn, method)
               }
               
               # "Mixed" gene-set statistics
               gsStatsAbs[iPerm] <- calcGeneSetStat(abs(randStats), method)
            }
            
            
            #************************
            # Wilcoxon:
            #************************
            if(method == "wilcoxon") {
             
               # Select random set:
               #rInd <- ceiling(runif(nGenesInSet,0,length(statsContrast)))
               rInd <- sample(1:length(statsContrast),nGenesInSet)
               randStats <- statsContrast[rInd]
               statsNotInSet <- statsContrast[!c(1:length(statsContrast))%in%rInd]
               
               # "Directional" gene-set statistics:
               if(statType == "p-signed") {
                  rInd <- sample(1:length(statsContrastTestUp),nGenesInSet)
                  randStatsTestUp <- statsContrastTestUp[rInd]
                  randStatsTestDn <- statsContrastTestDn[rInd]
                  statsNotInSetTestUp <- statsContrastTestUp[!c(1:length(statsContrastTestUp))%in%rInd]
                  statsNotInSetTestDn <- statsContrastTestDn[!c(1:length(statsContrastTestDn))%in%rInd]
                  gsStatsAllTestUp[iPerm] <- calcGeneSetStat(randStatsTestUp, "wilcoxon_fast", statsNotInSetTestUp)[1]
                  gsStatsAllTestDn[iPerm] <- calcGeneSetStat(randStatsTestDn, "wilcoxon_fast", statsNotInSetTestDn)[1]  
               } else if(statType == "t") {
                  gsStatsAll[iPerm] <- calcGeneSetStat(randStats, "wilcoxon_fast", statsNotInSet)[1]
               }
               
               # "Mixed" gene-set statistics: 
               gsStatsAbs[iPerm] <- calcGeneSetStat(abs(randStats), "wilcoxon_fast", abs(statsNotInSet))[1]
            }
            
            
            #************************
            # Mean:
            # Median:
            # Sum:
            #************************
            if(method %in% c("mean","median","sum")) {
               
               # Select random set:
               #rInd <- ceiling(runif(nGenesInSet,0,length(statsContrast)))
               rInd <- sample(1:length(statsContrast),nGenesInSet)
               randStats <- statsContrast[rInd]
               
               # "Directional" gene-set statistics:
               if(statType == "p-signed") {
                  rInd <- sample(1:length(statsContrastTestUp),nGenesInSet)
                  randStatsTestUp <- statsContrastTestUp[rInd]
                  randStatsTestDn <- statsContrastTestDn[rInd]
                  gsStatsAllTestUp[iPerm] <- calcGeneSetStat(randStatsTestUp, method)
                  gsStatsAllTestDn[iPerm] <- calcGeneSetStat(randStatsTestDn, method)
               } else if(statType == "t") {
                  gsStatsAll[iPerm] <- calcGeneSetStat(randStats, method)  
               }
               
               # "Mixed" gene-set statistics:
               gsStatsAbs[iPerm] <- calcGeneSetStat(abs(randStats), method)
               
            }
            
            
            #************************   
            # Maxmean:
            #************************
            if(method == "maxmean") {
               
               # "Mixed" gene-set statistics:
               #rInd <- ceiling(runif(nGenesInSet,0,length(statsContrast)))
               rInd <- sample(1:length(statsContrast),nGenesInSet)
               randStats <- statsContrast[rInd]
               gsStatsAbs[iPerm] <- calcGeneSetStat(randStats, method)
            }
               
            #************************   
            # GSEA:
            #************************
            if(method == "gsea") {
               
               # "Directional" gene-set statistics:
               #rInd <- ceiling(runif(nGenesInSet,0,length(statsContrast)))
               rInd <- sample(1:length(statsContrast),nGenesInSet)
               gsStatsAll[iPerm] <- calcGeneSetStat(rInd, method, statsContrastSorted, gseaParam)
            }
            
            
            #************************  
            # Page:
            #************************
            if(method == "page") {
             
               # "Directional" gene-set statistics:
               #rInd <- ceiling(runif(nGenesInSet,0,length(statsContrast)))
               rInd <- sample(1:length(statsContrast),nGenesInSet)
               randStats <- statsContrast[rInd]
               gsStatsAll[iPerm] <- calcGeneSetStat(randStats, method, statsContrast)
            }
           
         }
         
         # Save in matrix (each gene set size is represented by a row of statistics):
         gsStatsAllMatrix[iGeneSet,] <- gsStatsAll
         gsStatsAllTestUpMatrix[iGeneSet,] <- gsStatsAllTestUp
         gsStatsAllTestDnMatrix[iGeneSet,] <- gsStatsAllTestDn
         gsStatsAbsMatrix[iGeneSet,] <- gsStatsAbs
      }
      
      
      # For each gene set size (subset Up):
      if(statType %in% c("t","p-signed","F-signed") & method %in% c("fisher","stouffer","reporter","tailStrength","mean","median","sum","wilcoxon")) {
         for(iGeneSet in 1:length(gsSizesUp)) {
            gsStatsUp <- rep(NA,nPerm)
            nGenesInSet <- gsSizesUp[iGeneSet]
            
            # Loop nPerm times:
            if(nGenesInSet > 0) {
               for(iPerm in 1:nPerm) {
                  
                  #************************  
                  # Fisher:
                  # Stouffer:
                  # Reporter:
                  # Tail strength:
                  # Mean:
                  # Median:
                  # Sum:
                  #************************
                  if(method %in% c("fisher","stouffer","reporter","tailStrength","mean","median","sum")) {
                     #rInd <- ceiling(runif(nGenesInSet,0,length(statsContrastUp)))
                     rInd <- sample(1:length(statsContrastUp),nGenesInSet)
                     randStats <- statsContrastUp[rInd]
                     gsStatsUp[iPerm] <- calcGeneSetStat(randStats, method)
                     
                  } 
                  
                  #************************  
                  # Wilcoxon:
                  #************************
                  if(method == "wilcoxon") {
                     #rInd <- ceiling(runif(nGenesInSet,0,length(statsContrastUp)))
                     rInd <- sample(1:length(statsContrastUp),nGenesInSet)
                     randStats <- statsContrastUp[rInd]
                     statsNotInSet <- statsContrastUp[!c(1:length(statsContrastUp))%in%rInd]
                     gsStatsUp[iPerm] <- calcGeneSetStat(randStats, "wilcoxon_fast", statsNotInSet)[1]
                  }
               }
            }
            
            # Save in matrix (each gene set size is represented by a row of statistics):
            gsStatsUpMatrix[iGeneSet,] <- gsStatsUp
         }
      }
      
       
      # For each gene set size (subset Dn):
      if(statType %in% c("t","p-signed","F-signed") & method %in% c("fisher","stouffer","reporter","tailStrength","mean","median","sum","wilcoxon")) {
         for(iGeneSet in 1:length(gsSizesDn)) {
            gsStatsDn <- rep(NA,nPerm)
            nGenesInSet <- gsSizesDn[iGeneSet]
            
            # Loop nPerm times:
            if(nGenesInSet > 0) {
               for(iPerm in 1:nPerm) {
                  
                  #************************  
                  # Fisher:
                  # Stouffer:
                  # Reporter:
                  # Tail strength:
                  # Mean:
                  # Median:
                  # Sum:
                  #************************
                  if(method %in% c("fisher","stouffer","reporter","tailStrength","mean","median","sum")) {
                     #rInd <- ceiling(runif(nGenesInSet,0,length(statsContrastDn)))
                     rInd <- sample(1:length(statsContrastDn),nGenesInSet)
                     randStats <- statsContrastDn[rInd]
                     gsStatsDn[iPerm] <- calcGeneSetStat(randStats, method)
                     
                  } 
                  
                  #************************  
                  # Wilcoxon:
                  #************************
                  if(method == "wilcoxon") {
                     #rInd <- ceiling(runif(nGenesInSet,0,length(statsContrastDn)))
                     rInd <- sample(1:length(statsContrastDn),nGenesInSet)
                     randStats <- statsContrastDn[rInd]
                     statsNotInSet <- statsContrastDn[!c(1:length(statsContrastDn))%in%rInd]
                     gsStatsDn[iPerm] <- calcGeneSetStat(randStats, "wilcoxon_fast", statsNotInSet)[1]
                  }
               }
            }
            
            # Save in matrix (each gene set size is represented by a row of statistics):
            gsStatsDnMatrix[iGeneSet,] <- gsStatsDn
         }
      }
      
      # Save results:
      rownames(gsStatsAllMatrix) <- gsSizes
      rownames(gsStatsAllTestUpMatrix) <- gsSizes
      rownames(gsStatsAllTestDnMatrix) <- gsSizes
      rownames(gsStatsAbsMatrix) <- gsSizes
      rownames(gsStatsUpMatrix) <- gsSizesUp
      rownames(gsStatsDnMatrix) <- gsSizesDn
      res$gsStatsAllPerm[[iContrast]] <- gsStatsAllMatrix
      res$gsStatsAllTestUpPerm[[iContrast]] <- gsStatsAllTestUpMatrix
      res$gsStatsAllTestDnPerm[[iContrast]] <- gsStatsAllTestDnMatrix
      res$gsStatsAbsPerm[[iContrast]] <- gsStatsAbsMatrix
      res$gsStatsUpPerm[[iContrast]] <- gsStatsUpMatrix
      res$gsStatsDnPerm[[iContrast]] <- gsStatsDnMatrix
      
   }
   return(res)
   
}