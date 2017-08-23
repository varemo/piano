#' Heatmap of top consensus gene sets
#' 
#' Based on multiple result objects from the \code{\link{runGSA}} function,
#' this function computes the consensus scores, based on rank aggregation, for
#' each directionality class and produces a heatmap plot of the results.
#' 
#' This function computes the consensus gene set scores for each directionality
#' class based on the results (gene set p-values) listed in \code{resList},
#' using the \code{consensusScores} function. For each class, only the
#' \code{GSAres} objects in \code{resList} that contain p-values for that class
#' are used as a basis for the rank aggregation. Hence, if not all classes are
#' covered by at least 2 \code{GSAres} objects in the list, the
#' \code{consensusHeatmap} function will not work. The results are displayed in
#' a heatmap showing the consensus scores.
#' 
#' @param resList a list where each element is an object of class
#' \code{GSAres}, as returned by the \code{runGSA} function.
#' @param method a character string selecting the method, either \code{"mean"},
#' \code{"median"}, \code{"Borda"} or \code{"Copeland"}.
#' @param cutoff the maximum consensus score of a gene set, in any of the
#' directionality classes, to be included in the heatmap.
#' @param adjusted a logical, whether to use adjusted p-values or not. Note
#' that if \code{runGSA} was run with the argument \code{adjMethod="none"}, the
#' adjusted p-values will be equal to the original p-values.
#' @param plot whether or not to draw the heatmap. Setting \code{plot=FALSE}
#' allows you to save the heatmap as a matrix without plotting it.
#' @param ncharLabel the number of characters to include in the row labels.
#' @param cellnote a character string selecting the information to be printed
#' inside each cell of the heatmap. Either \code{"consensusScore"},
#' \code{"medianPvalue"}, \code{"nGenes"} or \code{"none"}. Note that the
#' actual heatmap will always be based on the consensus scores.
#' @param columnnames either \code{"full"} (default) or \code{"abbr"} to use
#' full or abbreviated column labels. Will save some space for the heatmap if
#' set to \code{"abbr"}
#' @param colorkey a logical (default \code{TRUE}), whether or not to display
#' the colorkey. Will save some space for the heatmap if turned off.
#' @param colorgrad a character vector giving the color names to use in the
#' heatmap.
#' @param cex a numeric, to control the text size.
#' @return A list, returned invisibly, containing the matrix of consensus
#' scores as represented in the heatmap as well as the matrix of corresponding
#' median p-values and the matrix of number of genes in each gene set
#' (inlcuding the subset of up and down regulated genes for the mixed
#' directional classes).
#' @author Leif Varemo \email{piano.rpkg@@gmail.com} and Intawat Nookaew
#' \email{piano.rpkg@@gmail.com}
#' @seealso \pkg{\link{piano}}, \code{\link{runGSA}}
#' @examples
#' 
#' 
#'    # Load some example GSA results:
#'    data(gsa_results)
#'    
#'    # Consensus heatmap:
#'    dev.new(width=10,height=10)
#'    consensusHeatmap(resList=gsa_results)
#'    
#'    # Store the output:
#'    dev.new(width=10,height=10)
#'    ch <- consensusHeatmap(resList=gsa_results)
#'    
#'    # Access the median p-values for gene set s1:
#'    ch$pMat["s1",]
#' 
consensusHeatmap <- function(resList, method="median", cutoff=5, adjusted=FALSE, plot=TRUE,  
                              ncharLabel=25, cellnote="consensusScore", columnnames="full",
                              colorkey=TRUE, colorgrad=NULL, cex=NULL) {
   
   # error check:
   tmp <- try(method <- match.arg(method, c("mean","median","Borda","Copeland"), several.ok=FALSE), silent=TRUE)
   if(class(tmp) == "try-error") {
      stop("argument method is not valid")
   }
   tmp <- try(cellnote <- match.arg(cellnote, c("consensusScore","medianPvalue","nGenes","none"), several.ok=FALSE), silent=TRUE)
   if(class(tmp) == "try-error") {
      stop("argument cellnote is not valid")
   }
   tmp <- try(columnnames <- match.arg(columnnames, c("full","abbr"), several.ok=FALSE), silent=TRUE)
   if(class(tmp) == "try-error") {
      stop("argument columnnames is not valid")
   }
   if(!is.null(colorgrad)) {
      if(!is.character(colorgrad) | length(colorgrad)<2) {
         stop("argument colorgrad should be a character vector of minimum length 2")
      }
   }
   if(length(cutoff) != 1 | cutoff < 1) stop("argument cutoff should be a positive integer")
   if(!is.logical(adjusted)) stop("argument adjusted should be a logical")
   
   # Get consensus ranks for each directionality class:
   reslist <- resList
   nod <- consensusScores(resList=reslist,class="non",adjusted=adjusted,method=method,
                          n=length(reslist[[1]]$gsc),plot=FALSE)
   ddu <- consensusScores(resList=reslist,class="distinct",adjusted=adjusted,direction="up",method=method,
                          n=length(reslist[[1]]$gsc),plot=FALSE)
   ddd <- consensusScores(resList=reslist,class="distinct",adjusted=adjusted,direction="down",method=method,
                          n=length(reslist[[1]]$gsc),plot=FALSE)
   mdu <- consensusScores(resList=reslist,class="mixed",adjusted=adjusted,direction="up",method=method,
                          n=length(reslist[[1]]$gsc),plot=FALSE)
   mdd <- consensusScores(resList=reslist,class="mixed",adjusted=adjusted,direction="down",method=method,
                          n=length(reslist[[1]]$gsc),plot=FALSE)
   
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
   ngenesmat <- cbind(nGenes,nGenesDn,nGenes,nGenesUp,nGenes)
   rownames(ngenesmat) <- rownames(plotmat)
   colnames(plotmat) <- c("Non-directional","Distinct-directional (up)","Distinct-directional (dn)",
                          "Mixed-directional (up)","Mixed-directional (dn)")
   #rownames(plotmat) <- paste(rownames(plotmat)," (Tot:",nGenes,", Up:",nGenesUp,", Down:",nGenesDn,")",sep="")
   myorder <- c(3,5,1,4,2)
   
   # median p-value matrix:
   gs_names <- rownames(plotmat[,myorder])
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
   
   if(plot) {
      set.seed(1) # <---------------------
      tmp <- round(rnorm(10000,0,1000)) # <------------------------- FIX COLORING BETTER!!!
      tmp <- rev(unique(sort(tmp[tmp>0])))
      if(is.null(colorgrad)) {
         clrs <- c("red3","red","orange","yellow","lightyellow","white")
      } else {
         clrs <- colorgrad  
      }
      mycol <- colorRampPalette(rev(clrs), interpolate = "linear")(max(tmp))[tmp]
      tmpMat <- plotmat
      tmp <- rownames(tmpMat)
      for(i in 1:length(tmp)) {
         if(nchar(tmp[i])>ncharLabel) tmp[i] <- paste(substr(tmp[i],1,ncharLabel),"...",sep="")
      }
      rownames(tmpMat) <- tmp
      labAngle <- 90
      labAdj <- NULL
      if(columnnames=="abbr") {
         colnames(tmpMat) <- c("Nondir","Dist(up)","Dist(dn)","Mix(up)","Mix(dn)")
         labAngle <- 0
         labAdj <- c(0.5,0)
      }
      
      if(is.null(cex)) {
         mycex <- 0.2 + 1/log10(nrow(tmpMat))
      } else {
         mycex <- cex  
      }
      
      if(cellnote=="consensusScore") {
         notemat <- round(tmpMat[,myorder],3)
      } else if(cellnote=="medianPvalue") {
         notemat <- pMat
      } else if(cellnote=="nGenes") {
         notemat <- ngenesmat  
      } else if(cellnote=="none") {
         notemat <- matrix(rep("",nrow(tmpMat)*ncol(tmpMat)),nrow(tmpMat),ncol(tmpMat))
      }
      topmar <- ifelse(colorkey,10,1)
      bottommar<- ifelse(columnnames=="full",10,0)
      
      #start new section
      #clustering by gene members
      #omat <- matrix(nrow=nrow(tmpMat),ncol=nrow(tmpMat))
      #for(i in 1:nrow(tmpMat)) {
      #   for(j in 1:nrow(tmpMat)) {
      #      tmp1 <- reslist[[1]]$gsc[[i]]
      #      tmp2 <- reslist[[1]]$gsc[[j]]
      #      omat[i,j] <- sum(tmp1%in%tmp2) / length(unique(c(tmp1,tmp2)))
      #   }
      #}
      #mydendr <- as.dendrogram(hclust(dist(omat)))
      #end new section
      
      hm2out <- heatmap.2(tmpMat[,myorder], Colv=FALSE, dendrogram="row", margins=c(1,1), density.info="none",
                          key=colorkey, trace="none", scale="none", cellnote=notemat, notecol="black",
                          col=mycol,notecex=mycex, srtCol=labAngle, adjCol=labAdj,
                          # change layout from 2*2 to 3*3:
                          lmat=matrix(c(3,2,7,4,1,8,5,6,9),nrow=3),
                          # set plot area heights according to number of rows:
                          lhei=c(topmar,   50,   25/nrow(tmpMat) + bottommar),
                          # set plot area width according to length of rownames and number of columns:
                          lwid=c(1,   ncol(tmpMat)/2,   0.8 + max(nchar(rownames(tmpMat)))/20),
                          # set col label text size to be the same as row text size:
                          cexCol=mycex, cexRow=mycex)
      hmRowInd <- rev(hm2out$rowInd)
   } else {
      hmRowInd <- 1:nrow(plotmat)
   }
   
   #Reorder pMat and ngenesmat to match heatmap
   pMat <- pMat[hmRowInd,]
   ngenesmat <- ngenesmat[hmRowInd,]
   
   # Return results:
   res <- list()
   res$rankMat <- plotmat[hmRowInd,myorder]
   res$pMat <- pMat
   res$nGenesMat <- ngenesmat
   invisible(res)
}










