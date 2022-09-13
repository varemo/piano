#' Gene set network plot
#' 
#' Draws a network with gene sets as nodes and the thickness of the edges
#' correlating to the number of shared genes. The gene set significance is
#' visualized as color intensities. Gives an overview of the influence of
#' overlap on significant gene sets.
#' 
#' In the case of \code{class="distinct"} and \code{direction="both"}, the
#' distinct directional p-values (\code{pDistinctDirUp} and
#' \code{pDistinctDirDn}, see \code{\link{runGSA}}) will be used in
#' combination. Using the \code{geneSets} and \code{lay} arguments, multiple
#' comparative plots (i.e. with the same layout) can be drawn, based for
#' instance on the output gene set list from other network plots with different
#' directionality classes.
#' 
#' @param gsaRes an object of class \code{GSAres}, as returned from
#' \code{runGSA()} or an object returned from \code{runGSAhyper()}.
#' @param class a character string determining the p-values of which
#' directionality class that should be used as significance information for the
#' plot. Can be one of \code{"distinct"}, \code{"mixed"}, \code{"non"}. Has to
#' be \code{"non"} if the result from \code{runGSAhyper()} is used.
#' @param direction a character string giving the direction of regulation, can
#' be either \code{"up"}, \code{"down"} or \code{"both"} (for
#' \code{class="distinct"} only).
#' @param adjusted a logical, if adjusted p-values should be used, or not. Note
#' that if \code{runGSA} was run with the argument \code{adjMethod="none"}, the
#' adjusted p-values will be equal to the original p-values.
#' @param significance the significance cut-off that determines which gene sets
#' are included in the plot. Defaults to 0.001.
#' @param geneSets a character vector of gene set names, to be included in the
#' plot. Defaults to \code{NULL}, but if given, the argument
#' \code{significance} will not be used.
#' @param overlap a positive numerical. Determines the smallest number of
#' sharing genes between two gene-sets that is needed in order to draw a
#' line/edge between the gene-sets. Defaults to 1.
#' @param lay a numerical between 1-5, or a layout function (see
#' \code{\link{layout}} in the \code{igraph} package). 1-5 sets the layout to
#' one of the five default layout for the network plot.
#' @param label a character string, either \code{"names"} ,\code{"numbers"},
#' \code{"numbersAndSizes"} or \code{"namesAndSizes"}, determining the labels
#' used for the nodes. The names are the gene set names, numbers is an
#' arbritary numbered list of the gene sets used in the plot connected to the
#' named list returned by the funtion. Sizes are the gene set sizes, e.g. the
#' number of genes.
#' @param cexLabel the text size of the node labels.
#' @param ncharLabel the number of characters to include in the node labels.
#' @param cexLegend the text size of the legend.
#' @param nodeSize a numerical vector of length 2 giving the maximum and
#' minimum node sizes. The node size represents the size of the gene set, and
#' all values will be scaled to the given interval.
#' @param edgeWidth a numerical vector of length 2 giving the maximum and
#' minimum edge widths. The edge width represents the number of shared genes
#' between two gene sets, and all values will be scaled to the given interval.
#' @param edgeColor a character vector giving the colors to use for increasing
#' edge width. Can also be set to a single color. Defaults to a gray-scale.
#' @param scoreColors a character vector giving the colors from which the
#' gradient used for node coloring will be created. In the case of
#' \code{class="distinct"} and \code{direction="both"} the first half of the
#' vector will be used for the up-regulated gene sets and the second part will
#' be used for the down-regulated gene sets.
#' @param main an optional character vector setting the title of the plot.
#' @return Returns a list with two components: \code{geneSets} containing the
#' names and numbers of the gene sets in the plot, and \code{layout},
#' containing the saved layout of the plot, which can be passed back to the
#' \code{lay} argument in order to draw a subsequent plot with the same layout.
#' @author Leif Varemo \email{piano.rpkg@@gmail.com} and Intawat Nookaew
#' \email{piano.rpkg@@gmail.com}
#' @seealso \pkg{\link{piano}}, \code{\link{runGSA}}, \code{\link{GSAheatmap}}, 
#' \code{\link{networkPlot2}}, \code{\link{exploreGSAres}},
#' \code{\link{layout}}
#' @examples
#' 
#'    # Load example input data to GSA:
#'    data("gsa_input")
#'    
#'    # Load gene set collection:
#'    gsc <- loadGSC(gsa_input$gsc)
#'       
#'    # Run gene set analysis:
#'    gsares <- runGSA(geneLevelStats=gsa_input$pvals , directions=gsa_input$directions, 
#'                     gsc=gsc, nPerm=500)
#'       
#'    # Network plot:
#'    networkPlot(gsares,class="non",significance=0.01)
#'    
#'    # Use circular layout and save the layout:
#'    nw <- networkPlot(gsares,class="non",significance=0.01,lay=5)
#'    
#'    # Use the saved layout to overlay the distinct-directional p-values for easy comparison.
#'    # Note that the gene sets are now not selected based on a significance cutoff, but from a list:
#'    networkPlot(gsares,class="distinct",direction="both",lay=nw$layout,geneSets=nw$geneSets)
#' 
networkPlot <- function(gsaRes, class, direction, adjusted=FALSE, significance=0.001, geneSets=NULL, overlap=1, lay=1, 
                        label="names", cexLabel=0.9, ncharLabel=25, cexLegend=1, nodeSize=c(10,40), edgeWidth=c(1,15), 
                        edgeColor=NULL, scoreColors=NULL, main) {
   
   message("-------------------------------------------------------------------------------\n
  A newer function 'networkPlot2()' with the same purpose is now available!\n
-------------------------------------------------------------------------------\n\n")
  
   test <- 1 # Which contrast? Currently only one allowed!
   
   #*********************************************
   # Check arguments:
   #*********************************************
   
   tmp <- try(pValue <- match.arg(class, c("distinct","mixed","non"), several.ok=FALSE), silent=TRUE)
   if(is(tmp, "try-error")) {
      stop("argument class is not valid")
   }
   if(pValue == "non") {
      if(!missing(direction)) warning("argument direction will not be used for pValue='non'")
      direction <- "none"
   } else {
      tmp <- try(direction <- match.arg(direction, c("up","down","both"), several.ok=FALSE), silent=TRUE)
      if(is(tmp, "try-error")) {
         stop("argument direction is not valid")
      }
   }
   
   if(pValue == "non") {
      maintext <- paste(pValue,"-directional",sep="")
   } else {
      maintext <- paste(pValue,"-directional (",direction,")",sep="")
   }
   
   
   if(pValue == "mixed" & direction == "both") stop("for pValue='mixed', direction can not be 'both'")
   
   if(pValue == "distinct" & direction == "both") pValue <- "dirupdn"
   if(pValue == "distinct" & direction == "up") pValue <- "dirup"
   if(pValue == "distinct" & direction == "down") pValue <- "dirdn"
   if(pValue == "non") pValue <- "mix"
   if(pValue == "mixed" & direction == "up") pValue <- "subup"
   if(pValue == "mixed" & direction == "down") pValue <- "subdn"
   
   #if(test < 1 | test > ncol(gsaRes$nGenes)) stop("argument test is out of bounds")
   if(significance < 0 | significance >1) stop("argument significance is out of bounds")
   if(overlap <= 0) stop("argument overlap has to be at least 1")
   if(length(nodeSize) != 2) stop("argument nodeSize has to have length 2")
   if(length(edgeWidth) != 2) stop("argument edgeWidth has to have length 2")
   if(!is(adjusted, "logical")) stop("argument adjusted has to be TRUE or FALSE")
   if(!missing(main)) if(!is(main, "character")) stop("argument main has to be a character string")
   
   #########################################################
   # Adds possibility to use output object from runGSAhyper:
   if(length(gsaRes) == 5) {
      if(all(names(gsaRes) == c("pvalues","p.adj","resTab","contingencyTable","gsc"))) {
         if(pValue != "mix") stop("When using result from runGSAhyper, only class='non' is allowed")
         gsaRes$pNonDirectional <- matrix(gsaRes$pvalues,ncol=1)
         gsaRes$pAdjNonDirectional <- matrix(gsaRes$p.adj,ncol=1) 
         gsaRes$geneSetStat <- "Fisher's exact test"
         if(adjusted) {
            maintext <- paste("p.adj<",significance,sep="")
         } else {
            maintext <- paste("p<",significance,sep="")
         }
      }
   }
   #########################################################
   
   
   #*********************************************
   # Prepare network:
   #*********************************************
   
   # Extract values:
   gsc <- gsaRes$gsc
   if(adjusted) { # ...if adjusted:
      if(pValue == "dirup") pValues <- gsaRes$pAdjDistinctDirUp[,test]
      if(pValue == "dirdn") pValues <- gsaRes$pAdjDistinctDirDn[,test]
      if(pValue == "dirupdn") {
         pValues <- apply(abs(cbind(gsaRes$pAdjDistinctDirUp[,test],gsaRes$pAdjDistinctDirDn[,test])),1,min,na.rm=TRUE)
         tmp     <- apply(abs(cbind(gsaRes$pAdjDistinctDirUp[,test],gsaRes$pAdjDistinctDirDn[,test])),1,which.min)==2
         pValues[pValues == 0] <- 1e-100
         pValues[tmp] <- -pValues[tmp]
      }
      if(pValue == "mix") pValues <- gsaRes$pAdjNonDirectional[,test]
      if(pValue == "subup") pValues <- gsaRes$pAdjMixedDirUp[,test]
      if(pValue == "subdn") pValues <- gsaRes$pAdjMixedDirDn[,test]
   
   } else { # ...if un-adjusted:
      if(pValue == "dirup") pValues <- gsaRes$pDistinctDirUp[,test]
      if(pValue == "dirdn") pValues <- gsaRes$pDistinctDirDn[,test]
      if(pValue == "dirupdn") {
         pValues <- apply(abs(cbind(gsaRes$pDistinctDirUp[,test],gsaRes$pDistinctDirDn[,test])),1,min,na.rm=TRUE)
         tmp     <- apply(abs(cbind(gsaRes$pDistinctDirUp[,test],gsaRes$pDistinctDirDn[,test])),1,which.min)==2
         pValues[pValues == 0] <- 1e-100
         pValues[tmp] <- -pValues[tmp]
      }
      if(pValue == "mix") pValues <- gsaRes$pNonDirectional[,test]
      if(pValue == "subup") pValues <- gsaRes$pMixedDirUp[,test]
      if(pValue == "subdn") pValues <- gsaRes$pMixedDirDn[,test]
   }
   
   # Get gene set names:
   geneSetNames <- names(gsc)
   
   # Alt. 1. Select user defined gene sets:
   if(!is.null(geneSets)) {
      indSignificant <- which(geneSetNames %in% geneSets)
      if(!missing(significance)) warning("argument significance will not be used when argument geneSets is supplied")
   
   # Alt. 2. Select significant gene sets:   
   } else {
      indSignificant <- which(abs(pValues) < significance)
   }
   
   # Check if at least two gene sets are selected:
   if(length(indSignificant) < 3) stop("less than three gene sets were selected, can not plot (tip: adjust the significance cutoff)")
   
   # Get the relevant p-values:
   pSignificant <- pValues[indSignificant]
   
   # Generate overlap matrix for significant gene sets:
   overlapMat <- matrix(nrow=length(indSignificant),ncol=length(indSignificant))
   for(i in 1:nrow(overlapMat)) {
      for(j in i:ncol(overlapMat)) {
         tmp <- sum(gsc[[indSignificant[i]]] %in% gsc[[indSignificant[j]]])
         overlapMat[i,j] <- tmp 
         overlapMat[j,i] <- tmp
      }
   }
   
   # Gene set size:
   gsSize <- diag(overlapMat)
   
   # Remove gene sets with only small overlap (defined by argument 'overlap'):
   overlapMat[overlapMat < overlap] <- 0
   
   # Adjecency matrix:
   adjMat <- overlapMat > 0
   
   # Create igraph object:
   tmp <- adjMat
   diag(tmp) <- 0 # For some reason diag=FALSE below stopped working...
   if(all(!tmp)) stop("no overlap between gene sets found, try to decrease argument overlap or increase argument significance")
   g <- graph.adjacency(tmp, mode="undirected", diag=FALSE)
   
   
   #*********************************************
   # Set plotting parameters:
   #*********************************************
   
   # Edge width, according to shared genes:
   edgeOverlap <- rep(NA,ecount(g))
   for(iEdge in 1:ecount(g)) {
      tmp <- ends(g,iEdge)
      edgeOverlap[iEdge] <- overlapMat[tmp[1],tmp[2]]
   }
   eWidth <- (edgeOverlap-min(edgeOverlap))/(max(edgeOverlap)-min(edgeOverlap))*(edgeWidth[2]-edgeWidth[1])+edgeWidth[1]
   
   # Edge color:
   if(is.null(edgeColor)) edgeColor = c("gray80","gray70","gray60")
   tmp <- seq(min(edgeOverlap), max(edgeOverlap), length.out=length(edgeColor)+1)
   eColor <- rep(edgeColor[1],ecount(g))
   for(i in 2:length(edgeColor)) {
      eColor[edgeOverlap > tmp[i]] <- edgeColor[i]  
   }
   
   # Node size, according to number of genes:
   vSize <- (gsSize-min(gsSize))/(max(gsSize)-min(gsSize))*(nodeSize[2]-nodeSize[1])+nodeSize[1] 
   
   # Node color:
   if(pValue == "dirupdn") {
      if(is.null(scoreColors)) {
         tmp1 <- c('red','tomato','mistyrose')
         tmp2 <- c('blue','cornflowerblue','azure')
         gradColorsUp <- colorRampPalette(tmp1,interpolate="linear")(100)
         gradColorsDn <- colorRampPalette(tmp2,interpolate="linear")(100)
      } else {
         if(length(scoreColors)%%2 != 0 | length(scoreColors) < 4) stop("argument scoreColors should contain at least four and an even number of colors")  
         tmp1 <- scoreColors[1:(length(scoreColors)/2)]
         tmp2 <- scoreColors[(length(scoreColors)/2+1):length(scoreColors)]
         gradColorsUp <- colorRampPalette(tmp1,interpolate="linear")(100)
         gradColorsDn <- colorRampPalette(tmp2,interpolate="linear")(100)
      }
      vColor <- rep(NA,length(pSignificant))
      tmp <- pSignificant[pSignificant > 0 & !is.na(pSignificant)]
      vColor[pSignificant > 0 & !is.na(pSignificant)] <- gradColorsUp[floor(tmp*99/max(tmp))+1]
      tmp <- abs(pSignificant[pSignificant < 0 & !is.na(pSignificant)])
      vColor[pSignificant < 0 & !is.na(pSignificant)] <- gradColorsDn[floor(tmp*99/max(tmp))+1]
   } else {
      if(is.null(scoreColors)) {
         tmp <- c('red','orange','yellow')
      } else {
         if(length(scoreColors) < 2) stop("argument scoreColors should contain at least two colors")
         tmp <- scoreColors
      }
      gradColors <- colorRampPalette(tmp,interpolate="linear")(100)
      if(max(pSignificant,na.rm=TRUE)==0) {
         vColor <- gradColors[rep(50,length(pSignificant))]   
      } else {
         vColor <- gradColors[floor(pSignificant*99/max(pSignificant,na.rm=TRUE))+1]  
      }
   } 
   vColor[is.na(vColor)] <- "#CCCCCC"
   
   # Node labels:
   tmp <- try(label <- match.arg(label, c("names","numbers","numbersAndSizes","namesAndSizes"), several.ok=FALSE), silent=TRUE)
   if(is(tmp, "try-error")) {
      stop("argument label has to be set to either 'names' or 'numbers'")
   }
   tmp <- names(gsc)[indSignificant]
   for(i in 1:length(tmp)) {
      if(nchar(tmp[i])>ncharLabel) tmp[i] <- paste(substr(tmp[i],1,ncharLabel),"...",sep="")
   }
   if(label == "names") vLabels <- tmp
   else if(label == "numbers") vLabels <- 1:length(indSignificant)
   else if(label == "numbersAndSizes") vLabels <- paste(1:length(indSignificant)," (",gsSize,")",sep="")
   else if(label == "namesAndSizes") vLabels <- paste(tmp," (",gsSize,")",sep="")
   
   
   #*********************************************
   # Set plotting layout:
   #*********************************************
   
   if(is(lay, "numeric")) {
      
      # Edge weight, inverse to number of edges and size (of the larger of the node pairs)
      # used by Fruchterman Reingold algorithm:
      eWeight <- rep(NA,ecount(g))
      for(iEdge in 1:ecount(g)) {
         tmp <- ends(g,iEdge)
         tmp1 <- gsSize[tmp[1]]
         tmp2 <- gsSize[tmp[2]]
         if(tmp1 > tmp2) {
            eWeight[iEdge] <- 1/(tmp1*sum(adjMat[,tmp[1]]))
         } else {
            eWeight[iEdge] <- 1/(tmp2*sum(adjMat[,tmp[2]]))
         }
      }
      if(min(eWeight) != max(eWeight)) {
         eWeight <- (eWeight-min(eWeight))/(max(eWeight)-min(eWeight))*(100-0)+0.1
      } else {
         eWeight <- rep(50,length(eWeight))
      }
      
      # Predefined layouts:
      if(lay == 1) lay <- layout_with_fr(g, weights=eWeight)      
      else if(lay == 2) lay <- layout_with_lgl(g)
      else if(lay == 3) lay <- layout_with_kk(g, weights=eWeight)
      else if(lay == 4) lay <- layout_with_graphopt(g, charge=vSize)
      else if(lay == 5) lay <- layout_in_circle(g)
   
   # User defined layout function:   
   } else if(is(lay, "function")){
      lay <- lay(g)  
   }
   

   #*********************************************
   # Plot:
   #*********************************************
   
   if(pValue == "dirupdn") {
      layout(matrix(c(1,4,2,4,3,3),ncol=3), widths=c(1,1,5), heights=c(2,3))
      par(cex.axis=cexLegend)
      maColorBar(seq(0,max(abs(pSignificant),na.rm=TRUE),by=max(abs(pSignificant),na.rm=TRUE)/100),FALSE,gradColorsUp, 
                 k=7, main="p-values up", cex.main=cexLegend)
      maColorBar(seq(0,max(abs(pSignificant),na.rm=TRUE),by=max(abs(pSignificant),na.rm=TRUE)/100),FALSE,gradColorsDn, 
                 k=7, main="p-values down", cex.main=cexLegend)
      
   } else {
      layout(matrix(c(1,3,2,2),ncol=2), widths=c(1,5), heights=c(1,1))
      maColorBar(seq(0,max(pSignificant,na.rm=TRUE),by=max(pSignificant,na.rm=TRUE)/100),FALSE,gradColors, 
                 k=7, main="p-values", cex.main=cexLegend)   
   }   
   if(missing(main)) main <- paste(gsaRes$geneSetStat,", ",maintext,sep="")
   plot(g, vertex.size=vSize, vertex.color=vColor, vertex.label=vLabels, vertex.label.family="sans",
        vertex.label.color="black", vertex.label.cex=cexLabel, edge.width=eWidth, 
        edge.color=eColor, layout=lay, main=main)
   plot.new()
   par(usr=c(-1,1,-1,1))
   text(0,0, cex=cexLegend,
        labels=paste("Max node size:\n",max(gsSize)," genes\n\nMin node size:\n",min(gsSize)," genes\n\n",
                     "Max edge width:\n",max(edgeOverlap)," genes\n\nMin edge width:\n",min(edgeOverlap),
                     " genes",sep=""), font=1, adj=0.5)
   
   # Reset layout:
   par(mfrow=c(1,1))
   
   
   
   #*********************************************
   # Return results:
   #*********************************************
   
   res <- list()
   res$geneSets <- names(gsc)[indSignificant]
   names(res$geneSets) <- 1:length(indSignificant)
   res$layout <- lay
   invisible(res)
   
}
