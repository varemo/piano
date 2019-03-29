networkPlot2 <- function(gsaRes, class, direction, adjusted=TRUE, significance=0.001, geneSets=NULL, overlap=0.1, 
                         lay="visNetwork", physics=T, label="names", labelSize=22, ncharLabel=25, nodeSize=c(10,40), 
                         edgeWidth=c(1,15), edgeColor=NULL, scoreColors=NULL, naColor="yellow", main, submain, seed=1,
                         maxAllowedNodes=Inf, shiny=FALSE) {
  
  test <- 1 # Which contrast? Currently only one allowed!
  
  #*********************************************
  # Check arguments:
  #*********************************************
  
  tmp <- try(pValue <- match.arg(class, c("distinct","mixed","non"), several.ok=FALSE), silent=TRUE)
  if(class(tmp) == "try-error") {
    stop("argument class is not valid")
  }
  if(pValue == "non") {
    if(!missing(direction)) warning("argument direction will not be used for pValue='non'")
    direction <- "none"
  } else {
    tmp <- try(direction <- match.arg(direction, c("up","down","both"), several.ok=FALSE), silent=TRUE)
    if(class(tmp) == "try-error") {
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
  if(overlap <= 0) stop("argument overlap has to be larger than zero")
  if(length(nodeSize) != 2) stop("argument nodeSize has to have length 2")
  if(length(edgeWidth) != 2) stop("argument edgeWidth has to have length 2")
  if(class(adjusted) != "logical") stop("argument adjusted has to be TRUE or FALSE")
  if(!missing(main)) if(!class(main) %in% c("character","NULL")) stop("argument main has to be a character string")
  
  if(!lay%in% c("visNetwork","layout_nicely","layout_as_star","layout_with_fr","layout_with_kk",
               "layout_with_sugiyama","layout_in_circle","layout_on_grid","layout_as_tree",
               "layout_on_sphere","layout_randomly","layout_with_dh","layout_with_gem",
               "layout_with_graphopt","layout_with_lgl","layout_with_mds",
               "1","2","3","4","5","6","7","8")) stop("layout not recognized")
  
  
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
      pValuesZero <- pValues == 0
      pValues[pValuesZero] <- min(c(pValues[pValues>0]/10, 1e-10), na.rm=TRUE)
      pValuesOne <- pValues == 1
      pValues[pValuesOne] <- 1/(1+1e-10)
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
      pValuesZero <- pValues == 0
      pValues[pValuesZero] <- min(c(pValues[pValues>0]/10, 1e-10), na.rm=TRUE)
      pValuesOne <- pValues == 1
      pValues[pValuesOne == 1] <- 1/(1+1e-10)
      pValues[tmp] <- -pValues[tmp]
    }
    if(pValue == "mix") pValues <- gsaRes$pNonDirectional[,test]
    if(pValue == "subup") pValues <- gsaRes$pMixedDirUp[,test]
    if(pValue == "subdn") pValues <- gsaRes$pMixedDirDn[,test]
  }
  if(pValue != "dirupdn") {
    pValuesZero <- pValues == 0
    pValues[pValuesZero] <- min(c(pValues[pValues>0]/10, 1e-10), na.rm=TRUE)
    pValuesOne <- rep(FALSE, length(pValues))
  }
  
  # Get gene set names:
  geneSetNames <- names(gsc)
  
  # Alt. 1. Select user defined gene sets:
  if(!is.null(geneSets)) {
    if(!all(geneSets %in% geneSetNames)) stop("argument geneSets not matching gene-set names in argument gsaRes")
    indSelected <- which(geneSetNames %in% geneSets)
    if(!missing(significance)) warning("argument significance will not be used when argument geneSets is supplied")
    
    # Alt. 2. Select significant gene sets:   
  } else {
    indSelected <- which(abs(pValues) < significance)
  }
  
  # Check if at least two gene sets are selected:
  if(length(indSelected) < 2) stop("less than two gene sets were selected, can not plot (tip: adjust the significance cutoff)")
  
  # Check if too many gene sets are selected:
  if(length(indSelected) > maxAllowedNodes) stop(paste("the selected parameters results in a network with more than ",maxAllowedNodes," nodes (gene-sets). Drawing large networks requires more memory, if you want to continue, adjust the maxAllowedNodes argument.",sep=""))
  
  # Get the relevant p-values and convert to -log10:
  pSelected <- pValues[indSelected]
  pSelectedLog10 <- -sign(pSelected)*log10(abs(pSelected))
  
  # Generate overlap matrix for significant gene sets:
  overlapMat <- matrix(nrow=length(indSelected),ncol=length(indSelected))
  for(i in 1:nrow(overlapMat)) {
    for(j in i:ncol(overlapMat)) {
      tmp <- sum(gsc[[indSelected[i]]] %in% gsc[[indSelected[j]]])
      overlapMat[i,j] <- tmp 
      overlapMat[j,i] <- tmp
    }
  }
  
  # Gene set size:
  gsSize <- diag(overlapMat)
  
  # Generate overlap percentage matrix (% of smallest geneset in pair) 
  overlapPercentMat <- matrix(nrow=length(indSelected),ncol=length(indSelected))
  for(i in 1:nrow(overlapPercentMat)) {
    for(j in i:ncol(overlapPercentMat)) {
      tmp <- overlapMat[i,j]/min(gsSize[i],gsSize[j])
      overlapPercentMat[i,j] <- tmp
      overlapPercentMat[j,i] <- tmp
    }
  }
  
  # Remove gene sets with only small overlap (defined by argument 'overlap'):
  if(overlap>=1) {
    overlapMat[overlapMat < overlap] <- 0
  } else {
    overlapMat[overlapPercentMat < overlap] <- 0
  }
  
  # Adjecency matrix:
  adjMat <- overlapMat > 0
  
  # Create igraph object:
  tmp <- adjMat
  diag(tmp) <- 0 # For some reason diag=FALSE below stopped working...
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
  if(is.null(edgeColor)) edgeColor = c("gray90","gray80","gray70","gray60","gray50","gray40")
  tmp <- seq(min(edgeOverlap), max(edgeOverlap), length.out=length(edgeColor)+1)
  eColor <- rep(edgeColor[1],ecount(g))
  for(i in 2:length(edgeColor)) {
    eColor[edgeOverlap > tmp[i]] <- edgeColor[i]
  }

  # Node size, according to number of genes:
  vSize <- (gsSize-min(gsSize))/(max(gsSize)-min(gsSize))*(nodeSize[2]-nodeSize[1])+nodeSize[1]
  
  # Node color:
  colorLegendInfo <- list()
  if(pValue == "dirupdn") {
    if(is.null(scoreColors)) {
      tmp1 <- c('white','mistyrose','tomato','red')
      tmp2 <- c('white','azure','cornflowerblue','blue')
      gradColorsUp <- colorRampPalette(tmp1,interpolate="linear")(100)
      gradColorsDn <- colorRampPalette(tmp2,interpolate="linear")(100)
    } else {
      if(length(scoreColors)%%2 != 0 | length(scoreColors) < 4) stop("argument scoreColors should contain at least four and an even number of colors")
      tmp1 <- scoreColors[1:(length(scoreColors)/2)]
      tmp2 <- scoreColors[(length(scoreColors)/2+1):length(scoreColors)]
      gradColorsUp <- colorRampPalette(tmp1,interpolate="linear")(100)
      gradColorsDn <- colorRampPalette(tmp2,interpolate="linear")(100)
    }
    vColor <- rep(NA,length(pSelectedLog10))
    tmp <- pSelectedLog10[pSelectedLog10 > 0 & !is.na(pSelectedLog10)]
    vColor[pSelectedLog10 > 0 & !is.na(pSelectedLog10)] <- gradColorsUp[round(rescale(tmp,c(1,100),c(-1e-8,max(abs(pSelectedLog10),3,na.rm=T))))]
    tmp <- abs(pSelectedLog10[pSelectedLog10 < 0 & !is.na(pSelectedLog10)])
    vColor[pSelectedLog10 < 0 & !is.na(pSelectedLog10)] <- gradColorsDn[round(rescale(tmp,c(1,100),c(-1e-8,max(abs(pSelectedLog10),3,na.rm=T))))]
    
    colorLegendInfo$colors <- c(rev(gradColorsDn),gradColorsUp)
    colorLegendInfo$range <- max(abs(pSelectedLog10),3,na.rm=T)*c(-1,1)
    
  } else {
    if(is.null(scoreColors)) {
      if(direction=="up") {
        tmp <- c('white','mistyrose','tomato','red')
      } else if(direction=="down") {
        tmp <- c('white','azure','cornflowerblue','blue')
      } else {
        tmp <- c('#DDFABA','#55D800')
      }
    } else {
      if(length(scoreColors) < 2) stop("argument scoreColors should contain at least two colors")
      tmp <- scoreColors
    }
    gradColors <- colorRampPalette(tmp,interpolate="linear")(100)
    vColor <- rep(NA, length(pSelectedLog10))
    vColor[!is.na(pSelectedLog10)] <- gradColors[round(rescale(pSelectedLog10[!is.na(pSelectedLog10)],c(1,100),c(-1e-8,max(pSelectedLog10,3,na.rm=T))))]
    
    colorLegendInfo$colors <- gradColors
    colorLegendInfo$range <- c(0,max(pSelectedLog10,3,na.rm=T))
  }
  vColor[is.na(vColor)] <- naColor

  # Node labels:
  tmp <- try(label <- match.arg(label, c("names","numbers","numbersAndSizes","namesAndSizes"), several.ok=FALSE), silent=TRUE)
  if(class(tmp) == "try-error") {
    stop("argument label has to be set to either 'names' or 'numbers'")
  }
  tmp <- names(gsc)[indSelected]
  for(i in 1:length(tmp)) {
    if(nchar(tmp[i])>ncharLabel) tmp[i] <- paste(substr(tmp[i],1,ncharLabel),"...",sep="")
  }
  if(label == "names") vLabels <- tmp
  else if(label == "numbers") vLabels <- 1:length(indSelected)
  else if(label == "numbersAndSizes") vLabels <- paste(1:length(indSelected)," (",gsSize,")",sep="")
  else if(label == "namesAndSizes") vLabels <- paste(tmp," (",gsSize,")",sep="")
  
  # Node hover text:
  if(shiny) {
    tmp <- paste("<a href='#'"," onclick='Shiny.onInputChange(",'"links_genesets_click", "',names(gsc)[indSelected],'"',");'",">",names(gsc)[indSelected],"</a>",sep="")
  } else {
    tmp <- names(gsc)[indSelected]
  }
  if(pValue=="dirupdn") {
    tmp2 <- paste("<br>Direction:",ifelse(sign(pSelected)<0,"Down","Up"))
  } else {
    tmp2 <- "" 
  }
  tmp3 <- pSelected
  tmp3[pValuesZero[indSelected]] <- 0
  tmp3[pValuesOne[indSelected]] <- 1*sign(tmp3[pValuesOne[indSelected]])
  vHoverText <- paste(tmp, "<br>p-value: ", format(abs(tmp3),scientific=T,digits=3),"<br>-log10(p-value): ",round(log10(abs(tmp3)),3),tmp2,"<br>Genes: ", gsSize, sep="")
  
  #*********************************************
  # Predefined layouts:
  #*********************************************
  
  if(missing(physics) & lay%in%c(2:8)) physics <- FALSE
  
  if(lay == 2) lay <- "layout_nicely"
  else if(lay == 3) lay <- "layout_as_star"
  else if(lay == 4) lay <- "layout_with_fr"
  else if(lay == 5) lay <- "layout_with_kk"
  else if(lay == 6) lay <- "layout_with_sugiyama"
  else if(lay == 7) lay <- "layout_in_circle"
  else if(lay == 8) lay <- "layout_on_grid"
  
  
  #*********************************************
  # Main title and subtitle:
  #*********************************************
 
  if(missing(main)) {
    if(gsaRes$geneSetStat == "Fisher's exact test") {
      main <- paste(gsaRes$geneSetStat,", ",maintext,sep="")
    } else {
      main <- paste("GSA method: ",gsaRes$geneSetStat,", p-value: ",maintext,sep="")
    }
  }
  if(!is.null(main)) {
    main <- list(text=main, style="font-family:Arial;font-weight:normal;font-size:20px;text-align:center;")
  }
  
  if(missing(submain)) {
    submain <- paste("Node size (min-max): ",min(gsSize),"-",max(gsSize)," genes. Edge width (min-max): ",min(edgeOverlap),"-",max(edgeOverlap)," genes.", sep="")
  }
  if(!is.null(submain)) {
    submain <- list(text=submain, style="font-family:Arial;font-weight:normal;font-size:14px;text-align:center;")  
  }
  
  
  #*********************************************
  # Interactive plotting with visNetwork:
  #*********************************************
  
  # Convert igraph object to visNetwork:
  vn <- toVisNetworkData(g)
  
  # Node attributes:
  vn$nodes$geneSetNames <- names(gsc)[indSelected]
  vn$nodes$size <- vSize
  vn$nodes$label <- vLabels
  vn$nodes$color <- vColor
  vn$nodes$title <- vHoverText
  
  # Edge attributes:
  if(nrow(vn$edges)>0) {
    vn$edges$width <- eWidth
    vn$edges$color <- unlist(lapply(eColor, function(x) { rgb(col2rgb(x)[1],col2rgb(x)[2],col2rgb(x)[3], maxColorValue=255)}))
    vn$edges$title <- paste("Genes:",edgeOverlap)
  }
  
  # Plot:
  if(lay=="visNetwork" | lay==1) {
    res <- visNetwork(nodes=vn$nodes, edges=vn$edges, main=main, submain=submain) %>%
      visNodes(font=list(size=as.character(labelSize),face="arial"), shadow=T) %>%
      visLayout(randomSeed=seed) %>%
      visPhysics(enabled=physics)
  } else {
    res <- visNetwork(nodes=vn$nodes, edges=vn$edges, main=main, submain=submain) %>%
      visNodes(font=list(size=as.character(labelSize),face="arial"), shadow=T) %>%
      visIgraphLayout(layout=lay, physics=physics, randomSeed=seed)
  }
  
  # Draw / optionally return object for later drawing:
  res$colorLegendInfo <- colorLegendInfo
  res
  
  # Examples of reusing res later:
  # Draw same again:
  # visNetwork(res$x$nodes,res$x$edges)
  # os simly just:
  # res
  # Draw only essential, rest is default:
  # visNetwork(res$x$nodes[,c("id","label")],res$x$edges[,c("from","to")])
  # Add custom options:
  # visNetwork(res$x$nodes[,c("id","label")],res$x$edges[,c("from","to")]) %>% visIgraphLayout("layout_in_circle")
  # Other example:
  # res %>% visNodes(shadow=F)
  # Display number to gene-set name mapping:
  # res$x$nodes[,c("id","geneSetNames")]
}



