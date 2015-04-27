polarPlot <- function(pValues, chromosomeMapping,
                      colors=c("red","green","blue","yellow","orange",
                               "purple","tan","cyan","gray60","black"),
                      save=FALSE, verbose=TRUE) {

  #if(!try(require(gtools))) stop("package gtools is missing") # old, line below is preferred:
  if (!requireNamespace("gtools", quietly = TRUE)) stop("package gtools is missing")
  #if(!try(require(plotrix))) stop("package plotrix is missing") # old, line below is preferred:
  if (!requireNamespace("plotrix", quietly = TRUE)) stop("package plotrix is missing")
  
  # Verbose function:
  .verb <- function(message, verbose) {
    if(verbose == TRUE) {
      message(message)
    }
  }
  
  savedirFig <- paste(getwd(),"/Piano_Results/Figures/DifferentialExpression", sep="")
   
  
  # Load chromosome mapping:
  if(class(chromosomeMapping) == "character") {
    if(file.exists(chromosomeMapping)) {
      chromosomeMapping <- as.data.frame(read.delim(chromosomeMapping, header=TRUE, sep="\t",
                           row.names=1, as.is=TRUE),stringsAsFactors=FALSE,quote="")
    } else {
      stop("could not find the annotation file")
    }
  } else {
    chromosomeMapping <- as.data.frame(chromosomeMapping)
  }
  
  
  nContrasts <- ncol(pValues)
  
  # error check
  if(nContrasts > 10) {
    stop("can only handle up to 10 comparisons")
  }
   
  chromosomeMapping[,1] <- as.character(chromosomeMapping[,1])
  # Remove minus sign from location (due to sense vs antisense)
  suppressWarnings(tmp <- as.numeric(as.character(chromosomeMapping[,2])))
  if(!all(!is.na(tmp))) stop("NAs in chromosomeMapping, chromosome location")
  chromosomeMapping[,2] <- abs(tmp)
   
   
  # Combine p-values and chr info
  allData <- merge(x=pValues,y=chromosomeMapping,by.x=0,by.y=0,all.x=TRUE,stringsAsFactors=FALSE)
  rownames(allData) <- allData[,1]
  allData <- allData[,c(2:ncol(allData))]
  colnames(allData)[c(nContrasts+1,ncol(allData))] <- c("chromosome","start")
  allData[is.na(allData[,nContrasts+1]),nContrasts+1] <- "U"
  allData[is.na(allData[,nContrasts+2]),nContrasts+2] <- 0
  
  chrNames <- c(gtools::mixedsort(unique(allData[,nContrasts+1])),"-log10")
  nChr <- length(chrNames)
   

  # Create plotData and radialPos with the line lengths and radial positions in plot order:
  plotData <- 0  # initialize
  radialPos <- 0  # initialize
  for(i in 1:nChr) {
    dataChrSubset <- allData[allData$chromosome == chrNames[i],]
    nGenesChr <- nrow(dataChrSubset)
    if(nGenesChr > 0) {
      ind <- order(dataChrSubset$start, na.last=TRUE)
      dataChrSubsetSorted <- dataChrSubset[ind,c(1:nContrasts)]

      if(nContrasts > 1) {
        plotData <- rbind(plotData,dataChrSubsetSorted)  # The lengths of the lines in plot order
      } else {
        plotData <- c(plotData,dataChrSubsetSorted)  # The lengths of the lines in plot order
      }

      radialPosChr <- seq(from=(i-1)*2*pi/nChr,to=i*2*pi/nChr,length.out=nGenesChr+1)
      radialPosChr <- radialPosChr[1:nGenesChr]
      radialPos <- c(radialPos,radialPosChr) # The radial positions of the lines in plot order
    }
  }
  radialPos <- radialPos[2:length(radialPos)]  # get rid of the first dummy-row
  if(nContrasts > 1) {
    plotData <- plotData[2:nrow(plotData),]  # get rid of the first dummy-row
  } else {
    plotData <- plotData[2:length(plotData)]
  }
  
  # Adjust the radial positions to a small rotation to allow for the axis
  radialPos <- radialPos + pi/nChr
  radialChrPos <- seq(from=0,to=2*pi,length.out=nChr+1) + pi/nChr
  radialChrLabPos <- radialChrPos + pi/nChr


  # Plot:
  if(save == TRUE) {
    dirStat <- dir.create(savedirFig, recursive=TRUE, showWarnings=FALSE)
    if(dirStat == TRUE) {
      .verb(paste("Creating new directory:",savedirFig), verbose)
    }
    .verb("Saving merged polar plot:...", verbose)
    pdfFilePath <- paste(savedirFig,"/","polarPlotMerged.pdf",sep="")
    if(file.exists(pdfFilePath)) {
      .verb("Warning: polarPlotMerged.pdf already exists in directory: overwriting old file...", verbose)
    }
  } else {
    .verb("Generating merged polar plot...", verbose)
  }
  maxRadius <- max(-log10(plotData))
  if(!is.vector(plotData)) {
    plotOrder <- sort(colMeans(plotData, na.rm=TRUE),index.retur=TRUE)$ix
  } else {
    plotOrder <- 1
  }
  # Add lines sectioning the circle into chromosome pie pieces
  if(save == TRUE) {
    pdf(file=pdfFilePath,paper="a4")
  } else {
    dev.new()
  }
  plotrix::radial.plot(c(0,rep(max(pretty(c(0,maxRadius))),nChr)*1.1),rp.type="r",
              radial.pos=radialChrPos,line.col="gray75", show.radial.grid=FALSE, 
              radlab=TRUE, labels=chrNames, label.pos=radialChrLabPos, clockwise=TRUE,
              show.grid=TRUE,radial.lim=c(0,maxRadius))
  # Add the first p-value lines
  if(nContrasts > 1) {
    lengths <- -log10(plotData[,plotOrder[1]])
  } else {
    lengths <- -log10(plotData)
  }
  plotrix::radial.plot(lengths, radial.pos=radialPos, rp.type="r", clockwise=TRUE, add=TRUE, 
              line.col=colors[plotOrder[1]])
  # Add additional p-values lines
  if(nContrasts > 1) {
    for(i in 2:ncol(plotData)) {
      lengths <- -log10(plotData[,plotOrder[i]])
      plotrix::radial.plot(lengths, radial.pos=radialPos, rp.type="r", clockwise=TRUE, add=TRUE, 
                  line.col=colors[plotOrder[i]])
    }
  }
  if(save == TRUE) {
    tmp <- dev.off()
  }
  .verb("...done", verbose)
  # Plot single separate figures:
  if(nContrasts > 1) {
    if(save == TRUE) {
      .verb("Saving separate polar plots...", verbose)
    } else {
      .verb("Generating separate polar plots...", verbose)
    }
    for(i in 1:nContrasts) {
      if(save == TRUE) {
        pdfFilePath <- paste(savedirFig,"/","polarPlot_",colnames(pValues)[i],".pdf",sep="")
        if(file.exists(pdfFilePath)) {
          .verb(paste("Warning: polarPlot_",i,".pdf already exists in directory: overwriting old file..."), verbose)
        }
      }
      if(save == TRUE) {
        pdf(file=pdfFilePath,paper="a4")
      } else {
        dev.new()
      }
      maxRadius <- max(-log10(plotData))
      plotrix::radial.plot(c(0,rep(max(pretty(c(0,maxRadius))),nChr)*1.1),rp.type="r",
                  radial.pos=radialChrPos,line.col="gray75", show.radial.grid=FALSE, 
                  radlab=TRUE, labels=chrNames, label.pos=radialChrLabPos, clockwise=TRUE,
                  show.grid=TRUE,radial.lim=c(0,maxRadius))
      # Add the first p-value lines
      lengths <- -log10(plotData[,plotOrder[i]])
      plotrix::radial.plot(lengths, radial.pos=radialPos, rp.type="r", clockwise=TRUE, 
                  add=TRUE, line.col=colors[i])
      if(save == TRUE) {
        tmp <- dev.off()
      }
    }
    .verb("...done", verbose)
  }
  # Plot legend:
  if(save == TRUE) {
    .verb("Saving polar plot legend...", verbose)
    pdfFilePath <- paste(savedirFig,"/","polarPlotLegend.pdf",sep="")
    if(file.exists(pdfFilePath)) {
      .verb("Warning: polarPlotLegend.pdf already exists in directory: overwriting old file...", verbose)
    }
    pdf(file=pdfFilePath,paper="a4")
  } else {
    .verb("Generating polar plot legend...", verbose)
    dev.new()
  }
  plot.new()
  title(main="Legend Polar plot")
  legend(x=0.05, y=0.8, legend=colnames(pValues),fill=colors)
  if(save == TRUE) {
    tmp <- dev.off()
  }
  .verb("...done", verbose)
}