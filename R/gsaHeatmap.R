gsaHeatmap <- function(gsaRes, cutoff=5, adjusted=FALSE, ncharLabel=25, cellnote="pvalue") {
   
   tmp <- try(cellnote <- match.arg(cellnote, c("pvalue","rank","nGenes","none"), several.ok=FALSE), silent=TRUE)
   if(class(tmp) == "try-error") {
      stop("argument cellnote is not valid")
   }
   if(cellnote=="pvalue") cellnote<-"medianPvalue"
   if(cellnote=="rank") cellnote<-"consensusScore"
   
   res <- consensusHeatmap(resList=list(gsaRes,gsaRes), method="median", cutoff=cutoff, adjusted=adjusted, 
                     ncharLabel=ncharLabel, plot=TRUE, cellnote=cellnote)
   invisible(res)
}