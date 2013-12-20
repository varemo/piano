GSAheatmap <- function(gsaRes, cutoff=5, adjusted=FALSE, ncharLabel=25, cellnote="pvalue", columnnames="full",
                       colorkey=TRUE, colorgrad=NULL, cex=NULL) {
   
   tmp <- try(cellnote <- match.arg(cellnote, c("pvalue","rank","nGenes","none"), several.ok=FALSE), silent=TRUE)
   if(class(tmp) == "try-error") {
      stop("argument cellnote is not valid")
   }
   if(cellnote=="pvalue") cellnote<-"medianPvalue"
   if(cellnote=="rank") cellnote<-"consensusScore"
   
   # This will allow the heatmap even with some columns NA. But need to fix this handling in 
   # consensusHeatmap and consensusScores first also....
   #if(all(is.na(gsaRes$pMixedDirUp[,1]))) gsaRes$pMixedDirUp[,1] <- rep(1,nrow(gsaRes$pMixedDirUp))
   #if(all(is.na(gsaRes$pMixedDirDn[,1]))) gsaRes$pMixedDirDn[,1] <- rep(1,nrow(gsaRes$pMixedDirDn))
   #if(all(is.na(gsaRes$pDistinctDirUp[,1]))) gsaRes$pDistinctDirUp[,1] <- rep(1,nrow(gsaRes$pDistinctDirUp))
   #if(all(is.na(gsaRes$pDistinctDirDn[,1]))) gsaRes$pDistinctDirDn[,1] <- rep(1,nrow(gsaRes$pDistinctDirDn))
   
   tmp <- FALSE
   if(all(is.na(gsaRes$pMixedDirUp[,1]))) tmp <- TRUE
   if(all(is.na(gsaRes$pMixedDirDn[,1]))) tmp <- TRUE
   if(all(is.na(gsaRes$pDistinctDirUp[,1]))) tmp <- TRUE
   if(all(is.na(gsaRes$pDistinctDirDn[,1]))) tmp <- TRUE
   if(tmp) stop("The heatmap will only be drawn if all p-value classes exist. This is not the case for the supplied gsaRes object. Take a look at the networkPlot function instead.")
   
   res <- consensusHeatmap(resList=list(gsaRes,gsaRes), method="median", cutoff=cutoff, adjusted=adjusted, 
                     ncharLabel=ncharLabel, plot=TRUE, cellnote=cellnote, columnnames=columnnames,
                           colorkey=colorkey, colorgrad=colorgrad, cex=cex)
   invisible(res)
}