consensusHeatmap <- function(resList, method="median", cutoff=5, adjusted=FALSE) {
   
   # error check:
   tmp <- try(method <- match.arg(method, c("mean","median","Borda","Copeland"), several.ok=FALSE), silent=TRUE)
   if(class(tmp) == "try-error") {
      stop("argument method is not valid")
   }
   if(length(cutoff) != 1 | cutoff < 1 | cutoff > length(resList[[1]]$gsc)) stop("argument cutoff should be a positive integer no larger than the number of gene sets")
   if(!is.logical(adjusted)) stop("argument adjusted should be a logical")
   
   # Get consensus ranks for each directionality class:
   reslist <- resList
   nod <- consensusScores(resList=reslist,class="non",method=method,n=length(reslist[[1]]$gsc),plot=FALSE)
   ddu <- consensusScores(resList=reslist,class="distinct",direction="up",method=method,n=length(reslist[[1]]$gsc),plot=FALSE)
   ddd <- consensusScores(resList=reslist,class="distinct",direction="down",method=method,n=length(reslist[[1]]$gsc),plot=FALSE)
   mdu <- consensusScores(resList=reslist,class="mixed",direction="up",method=method,n=length(reslist[[1]]$gsc),plot=FALSE)
   mdd <- consensusScores(resList=reslist,class="mixed",direction="down",method=method,n=length(reslist[[1]]$gsc),plot=FALSE)
   
   nodPval <- nod$pMat
   dduPval <- ddu$pMat
   dddPval <- ddd$pMat
   mduPval <- mdu$pMat
   mddPval <- mdd$pMat
   
   nod <- nod$rankMat
   ddu <- ddu$rankMat
   ddd <- ddd$rankMat
   mdu <- mdu$rankMat
   mdd <- mdd$rankMat
   
   nod <- cbind(rownames(nod),nod[,2])
   colnames(nod) <- c("","nod")
   ddu <- cbind(rownames(ddu),ddu[,2])
   colnames(ddu) <- c("","ddu")
   ddd <- cbind(rownames(ddd),ddd[,2])
   colnames(ddd) <- c("","ddd")
   mdu <- cbind(rownames(mdu),mdu[,2])
   colnames(mdu) <- c("","mdu")
   mdd <- cbind(rownames(mdd),mdd[,2])
   colnames(mdd) <- c("","mdd")
   
   tmp <- merge(nod,ddu,by=1)
   tmp <- merge(tmp,ddd,by=1)
   tmp <- merge(tmp,mdu,by=1)
   tmp <- merge(tmp,mdd,by=1)
   
   consRankMat <- tmp[,2:6] # This is consensus scores.
   rownames(consRankMat) <- tmp[,1]
   for(i in 1:ncol(consRankMat)) {
      consRankMat[,i] <- as.numeric(as.character(consRankMat[,i]))
   }
   consRankMat <- as.matrix(consRankMat)
   
   plotmat <- consRankMat[which(apply(consRankMat,1,min)<=cutoff),]
   
   i <- 1
   nGenes <- rep(NA,nrow(plotmat))
   nGenesUp <- rep(NA,nrow(plotmat))
   nGenesDn <- rep(NA,nrow(plotmat))
   for(met in rownames(plotmat)) {
      nGenes[i] <- reslist[[1]]$nGenesTot[names(reslist[[1]]$gsc)==met]
      nGenesUp[i] <- reslist[[1]]$nGenesUp[names(reslist[[1]]$gsc)==met]
      nGenesDn[i] <- reslist[[1]]$nGenesDn[names(reslist[[1]]$gsc)==met]
      i <- i+1
   }
   notemat <- cbind(nGenes,nGenes,nGenes,nGenesUp,nGenesDn)
   colnames(plotmat) <- c("Non-directional","Distinct-directional (up)","Distinct-directional (dn)",
                          "Mixed-directional (up)","Mixed-directional (dn)")
   #rownames(plotmat) <- paste(rownames(plotmat)," (Tot:",nGenes,", Up:",nGenesUp,", Down:",nGenesDn,")",sep="")
   myorder <- c(3,5,1,4,2)
   set.seed(1) # <---------------------
   tmp <- round(rnorm(10000,0,1000)) # <------------------------- FIX COLORING BETTER!!!
   tmp <- rev(unique(sort(tmp[tmp>0])))
   mycol <- colorRampPalette(rev(c("red3","red","orange","yellow","lightyellow","white")), interpolate = "linear")(max(tmp))[tmp]
   tmpMat <- plotmat
   tmp <- rownames(tmpMat)
   for(i in 1:length(tmp)) {
      if(nchar(tmp[i])>25) tmp[i] <- paste(substr(tmp[i],1,25),"...",sep="")
   }
   rownames(tmpMat) <- tmp
   hm2out <- heatmap.2(tmpMat[,myorder], Colv=FALSE, dendrogram="row", margins=c(1,1), density.info="none",
                       key=TRUE, trace="none", scale="none", cellnote=round(tmpMat[,myorder],3), notecol="black",
                       col=mycol,notecex=0.2 + 1/log10(nrow(tmpMat)),
                       # change layout from 2*2 to 3*3:
                       lmat=matrix(c(3,2,7,4,1,8,5,6,9),nrow=3),
                       # set plot area heights according to number of rows:
                       lhei=c(10,   50,   25/nrow(tmpMat) + 10),
                       # set plot area width according to length of rownames and number of columns:
                       lwid=c(1,   ncol(tmpMat)/2,   0.8 + max(nchar(rownames(tmpMat)))/20),
                       # set col label text size to be the same as row text size:
                       cexCol=0.2 + 1/log10(nrow(tmpMat)))
   
   
   # median p-value matrix:
   gs_names <- rownames(plotmat[,myorder])[rev(hm2out$rowInd)]
   pMat <- matrix(ncol=5,nrow=length(gs_names))
   colnames(pMat) <- c("Distinct-directional (dn)","Mixed-directional (dn)","Non-directional","Mixed-directional (up)",
                       "Distinct-directional (up)")
   rownames(pMat) <- gs_names
   
   iGs <- match(gs_names,rownames(dddPval))
   pMat[,1] <- apply(dddPval[iGs,],1,median,na.rm=TRUE)
   iGs <- match(gs_names,rownames(mddPval))
   pMat[,2] <- apply(mddPval[iGs,],1,median,na.rm=TRUE)
   iGs <- match(gs_names,rownames(nodPval))
   pMat[,3] <- apply(nodPval[iGs,],1,median,na.rm=TRUE)
   iGs <- match(gs_names,rownames(mduPval))
   pMat[,4] <- apply(mduPval[iGs,],1,median,na.rm=TRUE)
   iGs <- match(gs_names,rownames(dduPval))
   pMat[,5] <- apply(dduPval[iGs,],1,median,na.rm=TRUE)
   
   # Return results:
   res <- list()
   res$rankMat <- plotmat[rev(hm2out$rowInd),myorder]
   res$pMat <- pMat
   invisible(res)
}










