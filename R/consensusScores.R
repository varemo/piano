consensusScores <- function(resList, class, direction, n=50, adjusted=FALSE, method="median", plot=TRUE, 
                            cexLabel=0.8, cexLegend=1, showLegend=TRUE, rowNames="names", logScale=FALSE, main) {
   
   test <- 1 # Which contrast? Only one allowed, currently!
   
   # Error check:
   tmp <- try(pValue <- match.arg(class, c("distinct","mixed","non"), several.ok=FALSE), silent=TRUE)
   if(class(tmp) == "try-error") {
      stop("argument class is not valid")
   }
   if(pValue == "non") {
      if(!missing(direction)) warning("argument direction will not be used for pValue='non'")
      direction <- "none"
   } else {
      tmp <- try(direction <- match.arg(direction, c("up","down"), several.ok=FALSE), silent=TRUE)
      if(class(tmp) == "try-error") {
         stop("argument direction is not valid")
      }
   }
   if(pValue == "distinct" & direction == "up") pValue <- "dirup"
   if(pValue == "distinct" & direction == "down") pValue <-"dirdn"
   if(pValue == "non") pValue <- "mix"
   if(pValue == "mixed" & direction == "up") pValue <-"subup"
   if(pValue == "mixed" & direction == "down") pValue <-"subdn"
   
   tmp <- try(method <- match.arg(method, c("mean","median","Borda","Copeland"), several.ok=FALSE), silent=TRUE)
   if(class(tmp) == "try-error") {
      stop("argument method is not valid")
   }
   tmp <- try(rowNames <- match.arg(rowNames, c("names","ranks","none"), several.ok=FALSE), silent=TRUE)
   if(class(tmp) == "try-error") {
      stop("argument rowNames is not valid")
   }
   
   # Check reslist length:
   if(length(resList)<2) stop("argument reslist contains too few objects")
   
   # Check similarity of results:
   tmp1 <- names(resList[[1]]$gsc)
   for(i in 2:length(resList)) {
      tmp2 <- names(resList[[i]]$gsc)
      if(length(tmp1) != length(tmp2)) stop("the listed GSAres objects have to have the same number of gene-sets in the same order")
      if(!all(tmp1 == tmp2)) stop("the listed GSAres objects have to have the same number of gene-sets in the same order")
      tmp1 <- tmp2
   }
   # Check that results have the selected p-value:
   saveInd <- rep(TRUE, length(resList))
   for(i in 1:length(resList)) {
      a <- resList[[i]]
      #if(test < 1 | test > ncol(a$nGenes)) stop("argument test is out of bounds")
      if(!adjusted) {
         if(pValue == "dirup") if(all(is.na(resList[[i]]$pDistinctDirUp[,test]))) saveInd[i] <- FALSE
         if(pValue == "dirdn") if(all(is.na(resList[[i]]$pDistinctDirDn[,test]))) saveInd[i] <- FALSE
         if(pValue == "mix")   if(all(is.na(resList[[i]]$pNonDirectional[,test]))) saveInd[i] <- FALSE
         if(pValue == "subup") if(all(is.na(resList[[i]]$pMixedDirUp[,test]))) saveInd[i] <- FALSE
         if(pValue == "subdn") if(all(is.na(resList[[i]]$pMixedDirDn[,test]))) saveInd[i] <- FALSE
      } else {
         if(pValue == "dirup") if(all(is.na(resList[[i]]$pAdjDistinctDirUp[,test]))) saveInd[i] <- FALSE
         if(pValue == "dirdn") if(all(is.na(resList[[i]]$pAdjDistinctDirDn[,test]))) saveInd[i] <- FALSE
         if(pValue == "mix")   if(all(is.na(resList[[i]]$pAdjNonDirectional[,test]))) saveInd[i] <- FALSE
         if(pValue == "subup") if(all(is.na(resList[[i]]$pAdjMixedDirUp[,test]))) saveInd[i] <- FALSE
         if(pValue == "subdn") if(all(is.na(resList[[i]]$pAdjMixedDirDn[,test]))) saveInd[i] <- FALSE
      }
   }
   if(sum(saveInd) < 2) stop("Less than two objects in resList contain the selected p-value")
   resList <- resList[saveInd]
   if(sum(!saveInd)) {
      warning(paste("No values for selected p-value for resList[c(",paste(which(!saveInd),collapse=","),
                 ")], omitting those objects from the list",sep=""))
   }
   
   # Rank matrix:
   rankMat <- matrix(nrow=length(resList),ncol=length(names(resList[[1]]$gsc)))
   for(i in 1:length(resList)) {
      a <- resList[[i]]
      #if(test < 1 | test > ncol(a$nGenesTot)) stop("argument test is out of bounds")
      if(!adjusted) {
         if(pValue == "dirup") a <- a$pDistinctDirUp[,test]
         if(pValue == "dirdn") a <- a$pDistinctDirDn[,test]
         if(pValue == "mix") a <- a$pNonDirectional[,test]
         if(pValue == "subup") a <- a$pMixedDirUp[,test]
         if(pValue == "subdn") a <- a$pMixedDirDn[,test]
      } else {
         if(pValue == "dirup") a <- a$pAdjDistinctDirUp[,test]
         if(pValue == "dirdn") a <- a$pAdjDistinctDirDn[,test]
         if(pValue == "mix") a <- a$pAdjNonDirectional[,test]
         if(pValue == "subup") a <- a$pAdjMixedDirUp[,test]
         if(pValue == "subdn") a <- a$pAdjMixedDirDn[,test]
      }
      a[is.na(a)] <- 1                         # <------ good to do this? Yep, otherwise NA will get increasing ranks.
      rankMat[i,] <- rank(a,ties.method="min") # <------ min or average? Probably min is more fair.
   }
   
   # Consensus score for each gene-set:
   rankScore <- rep(NA,ncol(rankMat))
   if(method == "median") {
      for(iGS in 1:ncol(rankMat)) {
         rankScore[iGS] <- median(rankMat[,iGS])
      }
   } else if(method == "mean") {
      for(iGS in 1:ncol(rankMat)) {
         rankScore[iGS] <- mean(rankMat[,iGS])
      }
   } else {
      tmp <- as.data.frame(t(rankMat))
      rownames(tmp) <- 1001:(nrow(tmp)+1000)
      rankScore <- relation_class_ids(relation_consensus(tmp,method))
      names(rankScore) <- NULL
   }
   
   # Get top ranked gene-sets:
   topInd <- which(rank(rankScore,ties.method="min") <= n) # index of top gene sets (no order)
   topInd <- topInd[sort(rankScore[topInd],index.return=TRUE)$ix] # index in rank order
   topInd <- rev(topInd)
   
   # Get consensus ranks:
   consRanks <- rank(rankScore[topInd],ties.method="min")
   
   # Row names:
   topGSnames <- names(resList[[1]]$gsc)[topInd]
   if(rowNames == "names") {
      tmp <- topGSnames
      for(i in 1:length(tmp)) {
         if(nchar(tmp[i])>35) tmp[i] <- paste(substr(tmp[i],1,35),"...",sep="")
      }
      rowNames <- tmp
      yaxt <- "s"
      ylabel <- ""
   } else if(rowNames == "ranks") {
      rowNames <- consRanks
      yaxt <- "s"
      ylabel <- "Consensus ranks"
   } else {
      rowNames <- rep("",length(topInd))  
      yaxt <- "n"
      ylabel <- ""
   }
   
   # Plot:
   if(plot) {
      logScale <- ifelse(logScale,"x","") 
      if(missing(main)) main <- paste(method, " rank, based on p",pValue,sep="")
      if(yaxt == "n") {
         layout(matrix(c(1,2),ncol=2),widths=c(4,1))
         cexPoint <- 0.9
      } else {
         layout(matrix(c(0,1,2),ncol=3),widths=c(1,4,1))
         cexPoint <- 1.2
      }
      boxplot(rankMat[,topInd], range=1.5, log=logScale, col="black",border="red",outline=!showLegend,
              pars=list(whisklty=0,staplelty=0,boxlty=0,outpch=20,outcol="darkgray",outcex=0.4,medlwd=1.5),horizontal=T,names=rowNames,
              oma=c(0,0,0,20),mar=c(0,0,0,20), las=1, cex.axis=cexLabel, main=main,yaxt=yaxt,xlab="Individual ranks",
              ylab=ylabel,cex.lab=cexLabel)
      if(showLegend) {
         for(i in 1:ncol(rankMat[,topInd])) {
            points(rankMat[,topInd][,i],rep(i,nrow(rankMat)),col=rainbow(nrow(rankMat)),pch=15:18,lwd=1,cex=cexPoint)     
         }
      }
      boxplot(rankMat[,topInd],col="black", log=logScale, range=1.5,border="red",outline=FALSE,
              pars=list(whisklty=1,whisklwd=0.5,whiskcol="black",staplelty=0,boxlty=0,medlwd=1.5),
              add=TRUE,horizontal=T,names=rowNames,las=1, cex.axis=cexLabel,yaxt=yaxt,xlab="Individual ranks",
              ylab=ylabel,cex.lab=cexLabel)
   }
   if(is.null(names(resList))) {
      methodNames <- lapply(resList, function(x) x$geneSetStat)
   } else {
      methodNames <- names(resList)
   }
   if(plot & showLegend) {
      par(xpd=NA)
      x <- ifelse(logScale=="x",10^par("usr")[2]*1.2,par("usr")[2]*1.05)
      legend(x=x,y=par("usr")[4],legend=methodNames,pch=15:18,col=rainbow(nrow(rankMat)),cex=cexLegend)
   }
   
   # Individual p-values for the gene sets:
   gsNames <- rev(topGSnames)
   iGs <- match(gsNames,names(resList[[1]]$gsc))
   pMat <- matrix(nrow=length(iGs),ncol=length(resList))
   rownames(pMat) <- gsNames
   colnames(pMat) <- methodNames
   for(i in 1:length(resList)) {
      if(!adjusted) {
         if(pValue == "dirup") pMat[,i] <- resList[[i]]$pDistinctDirUp[iGs,test]
         if(pValue == "dirdn") pMat[,i] <- resList[[i]]$pDistinctDirDn[iGs,test]
         if(pValue == "mix")   pMat[,i] <- resList[[i]]$pNonDirectional[iGs,test]
         if(pValue == "subup") pMat[,i] <- resList[[i]]$pMixedDirUp[iGs,test]
         if(pValue == "subdn") pMat[,i] <- resList[[i]]$pMixedDirDn[iGs,test]
      } else {
         if(pValue == "dirup") pMat[,i] <- resList[[i]]$pAdjDistinctDirUp[iGs,test]
         if(pValue == "dirdn") pMat[,i] <- resList[[i]]$pAdjDistinctDirDn[iGs,test]
         if(pValue == "mix")   pMat[,i] <- resList[[i]]$pAdjNonDirectional[iGs,test]
         if(pValue == "subup") pMat[,i] <- resList[[i]]$pAdjMixedDirUp[iGs,test]
         if(pValue == "subdn") pMat[,i] <- resList[[i]]$pAdjMixedDirDn[iGs,test]
      }
   }
   
   # Return result:
   rankMat <- t(rankMat[,rev(topInd)])
   rownames(rankMat) <- rev(topGSnames)
   colnames(rankMat) <- methodNames
   rankMat <- cbind(rev(consRanks),rev(rankScore[topInd]),rankMat)
   colnames(rankMat)[1:2] <- c("ConsRank","ConsScore")
   res <- list()
   res$rankMat <- rankMat
   res$pMat <- pMat
   return(res)
   
}
