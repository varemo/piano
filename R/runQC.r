runQC <- function(arrayData, rnaDeg=TRUE, nuseRle=TRUE, hist=TRUE,
                           boxplot=TRUE, pca=TRUE, colorFactor=1,
                           colors=c("red","green","blue","yellow","orange",
                                    "purple","tan","cyan","gray60","black",
                                    "white"),
                           save=FALSE, verbose=TRUE) {

   #if(!try(require(affyPLM))) stop("package affyPLM is missing") # old, line below is preferred:
   if (!requireNamespace("affyPLM", quietly = TRUE)) stop("package affyPLM is missing")
   #if(!try(require(affy))) stop("package affy is missing") # old, line below is preferred:
   if (!requireNamespace("affy", quietly = TRUE)) stop("package affy is missing")
  

  savedir <- paste(getwd(),"/Piano_Results/Figures/QualityControl", sep="")

  # First the functions for the various QCs (the code that calls the functions
  # is below):

  # Function for rna degradation:
  .runRnaDeg <- function(arrayData, savedir, save) {
    .verb("Generating RNA degradation plot...", verbose)
    rna.deg.obj <- affy::AffyRNAdeg(arrayData$dataRaw,log.it=TRUE)
    # Plot file:
    if(save == FALSE) {
      dev.new()
      affy::plotAffyRNAdeg(rna.deg.obj,transform="shift.scale")
    }
    .verb("...done", verbose)
    # Save file:
    if(save == TRUE) {
      dirStat <- dir.create(savedir, recursive=TRUE, showWarnings=FALSE)
      if(dirStat == TRUE) {
        .verb(paste("Creating new directory:",savedir), verbose)
      }
      .verb("Saving RNA degradation plot...", verbose)
      pdfFilePath = paste(savedir,"/","rnaDeg.pdf",sep="")
      if(file.exists(pdfFilePath)) {
        .verb("Warning: rnaDeg.pdf already exists in directory: overwriting old file...", verbose)
      }
      pdf(file=pdfFilePath,paper="a4")
      affy::plotAffyRNAdeg(rna.deg.obj,transform="shift.scale")
      tmp <- dev.off()
      .verb("...done", verbose)
    }
  }
  
  
  # Function for affyPLM:
  .runNuseRle <- function(arrayData, savedir, save) {
    .verb("Generating NUSE and RLE plot...", verbose)
    Pset <- affyPLM::fitPLM(arrayData$dataRaw)
    # Plot file:
    if(save == FALSE) {
      dev.new()
      par(mar=c(10, 4, 4, 2))
      #Mbox(Pset, main="RLE (Relative log expression)", las=2)
      affyPLM::RLE(Pset, main="RLE (Relative log expression)", las=2)
      dev.new()
      par(mar=c(10, 4, 4, 2))
      #boxplot(Pset, main="NUSE (Normalized unscaled standard error)", las=2)
      affyPLM::NUSE(Pset, main="NUSE (Normalized unscaled standard error)", las=2)
    }
    .verb("...done", verbose)
    # Save file:
    if(save == TRUE) {
      dirStat <- dir.create(savedir, recursive=TRUE, showWarnings=FALSE)
      if(dirStat == TRUE) {
        .verb(paste("Creating new directory:",savedir), verbose)
      }

      .verb("Saving RLE plot...", verbose)
      pdfFilePath = paste(savedir,"/","rle.pdf",sep="")
      if(file.exists(pdfFilePath)) {
        .verb("Warning: rle.pdf already exists in directory: overwriting old file...", verbose)
      }
      pdf(file=pdfFilePath,paper="a4")
      par(mar=c(10, 4, 4, 2))
      #Mbox(Pset, main="RLE (Relative Log Expression)", las=2)
      affyPLM::RLE(Pset, main="RLE (Relative Log Expression)", las=2)
      tmp <- dev.off()
      .verb("...done", verbose)
      
      .verb("Saving NUSE plot...", verbose)
      pdfFilePath = paste(savedir,"/","nuse.pdf",sep="")
      if(file.exists(pdfFilePath)) {
        .verb("Warning: nuse.pdf already exists in directory: overwriting old file...", verbose)
      }
      pdf(file=pdfFilePath,paper="a4")
      par(mar=c(10, 4, 4, 2))
      #boxplot(Pset, main="NUSE (Normalized unscaled standard error)", las=2)
      affyPLM::NUSE(Pset, main="NUSE (Normalized unscaled standard error)", las=2)
      tmp <- dev.off()
      .verb("...done", verbose)
    }
    
  }


  # Function for distribution:
  .runHist <- function(arrayData, savedir, save) {
    # Plot file:
    if(save == FALSE) {
      if("dataRaw" %in% attributes(arrayData)$names) {
        .verb("Generating raw data distribution plot...", verbose)
        dev.new()
        hist(arrayData$dataRaw, main="Histogram of raw data",  col="black",
             xlab=expression(log[2] ~ intensity), ylab="Density", lty=1)
        .verb("...done", verbose)
      }
      .verb("Generating normalized data distribution plot...", verbose)
      dev.new()
      plot(density(arrayData$dataNorm[,1]), main="Histogram of normalized data",
           xlab=expression(log[2] ~ intensity),ylab="Density")
      for(i in 2:dim(arrayData$dataNorm)[2]){
        lines(density(arrayData$dataNorm[,i]))
      }
      .verb("...done", verbose)
    }
    # Save file:
    if(save == TRUE) {
      dirStat <- dir.create(savedir, recursive=TRUE, showWarnings=FALSE)
      if(dirStat == TRUE) {
        .verb(paste("Creating new directory:",savedir), verbose)
      }

      if("dataRaw" %in% attributes(arrayData)$names) {
        .verb("Saving raw data distribution plot...", verbose)
        pdfFilePath = paste(savedir,"/","rawDataDist.pdf",sep="")
        if(file.exists(pdfFilePath)) {
          .verb("Warning: rawDataDist.pdf already exists in directory: overwriting old file...", verbose)
        }
        pdf(file=pdfFilePath,paper="a4")
        hist(arrayData$dataRaw, main="Histogram of raw data",  col="black",
             xlab=expression(log[2] ~ intensity), ylab="Density", lty=1)
        tmp <- dev.off()
        .verb("...done", verbose)
      }

      .verb("Saving normalized data distribution plot...", verbose)
      pdfFilePath = paste(savedir,"/","normDataDist.pdf",sep="")
      if(file.exists(pdfFilePath)) {
        .verb("Warning: normDataDist.pdf already exists in directory: overwriting old file...", verbose)
      }
      pdf(file=pdfFilePath,paper="a4")
      plot(density(arrayData$dataNorm[,1]), main="Histogram of normalized data",
           xlab=expression(log[2] ~ intensity),ylab="Density")
      for(i in 2:dim(arrayData$dataNorm)[2]){
        lines(density(arrayData$dataNorm[,i]))
      }
      tmp <- dev.off()
      .verb("...done", verbose)
    }
  }


  # Function for boxplot:
  .runBoxplot <- function(arrayData, savedir, save) {
    # Plot file:
    if(save == FALSE) {
      if("dataRaw" %in% attributes(arrayData)$names) {
        .verb("Generating raw data boxplot...", verbose)
        dev.new()
        par(mar=c(10, 4, 4, 2))
        boxplot(arrayData$dataRaw, main="Boxplot of raw data",
                ylab=expression(log[2] ~ intensity), las=2)
        .verb("...done", verbose)
      }
      .verb("Generating normalized data boxplot...", verbose)
      dev.new()
      par(mar=c(10, 4, 4, 2))
      boxplot(arrayData$dataNorm, main="Boxplot of normalized data",
                ylab=expression(log[2] ~ intensity), las=2)
      .verb("...done", verbose)
    }
    # Save file:
    if(save == TRUE) {
      dirStat <- dir.create(savedir, recursive=TRUE, showWarnings=FALSE)
      if(dirStat == TRUE) {
        .verb(paste("Creating new directory:",savedir), verbose)
      }

      if("dataRaw" %in% attributes(arrayData)$names) {
        .verb("Saving raw data boxplot...", verbose)
        pdfFilePath = paste(savedir,"/","rawDataBoxplot.pdf",sep="")
        if(file.exists(pdfFilePath)) {
          .verb("Warning: rawDataBoxplot.pdf already exists in directory: overwriting old file...", verbose)
        }
        pdf(file=pdfFilePath,paper="a4")
        par(mar=c(10, 4, 4, 2))
        boxplot(arrayData$dataRaw, main="Boxplot of raw data",
                ylab=expression(log[2] ~ intensity), las=2)
        tmp <- dev.off()
        .verb("...done", verbose)
      }

      .verb("Saving normalized data boxplot...", verbose)
      pdfFilePath = paste(savedir,"/","normDataBoxplot.pdf",sep="")
      if(file.exists(pdfFilePath)) {
        .verb("Warning: normDataBoxplot.pdf already exists in directory: overwriting old file...", verbose)
      }
      pdf(file=pdfFilePath,paper="a4")
      par(mar=c(10, 4, 4, 2))
      boxplot(arrayData$dataNorm, main="Boxplot of normalized data",
                ylab=expression(log[2] ~ intensity), las=2)
      tmp <- dev.off()
      .verb("...done", verbose)
    }
  }


  # Function for PCA:
  .runPca <- function(arrayData, savedir, save, colorFactor, colors) {
    .verb("Generating PCA...", verbose)
    dataForPca <- arrayData$dataNorm
    # Centralize data
    dataForPca <- dataForPca - rowMeans(dataForPca)
    
    # PCA
    dataPrcomp <- prcomp(dataForPca)
    # Order by pc3 (to plot in size order)
    pc3 <- dataPrcomp$rotation[,3]
    sizeOrder <- order(pc3, decreasing=TRUE)
    pc3 <- pc3[sizeOrder]
    pc1 <- dataPrcomp$rotation[,1]
    pc1 <- pc1[sizeOrder]
    pc2 <- dataPrcomp$rotation[,2]
    pc2 <- pc2[sizeOrder]
    # Scale pc3 size to better interval
    MinSize <- 1.5
    MaxSize <- 4.5
    pc3Size <- MinSize + (MaxSize - MinSize) * (pc3 - min(pc3))/(max(pc3) - min(pc3))
    # Coloring
    Colors <- colors
    mainFactors <- unique(arrayData$setup[,colorFactor])
    colorKey <- arrayData$setup[,colorFactor]
    colorKeyFactors <- mainFactors
    for(i in 1:length(mainFactors)) {
      colorKey[colorKey == mainFactors[i]] <- Colors[i]
      colorKeyFactors[colorKeyFactors == mainFactors[i]] <- Colors[i]
    }
    orderIndex <- NA
    for(i in 1:length(names(pc1))) {
      orderIndex[i] <- which(rownames(arrayData$setup) == names(pc1)[i])
    }
    colorKey <- colorKey[orderIndex]
    .verb("...done", verbose)

    # Plot:
    if(save == FALSE) {
      # Variance
      dev.new()
      mp <- barplot(summary(dataPrcomp)$importance[2,],names.arg=colnames(dataPrcomp$rotation),
                    las=3,ylab="Proportion of variance",
                    ylim=c(0,1.25*summary(dataPrcomp)$importance[2,][1]),
                    main="PC importance")
      text(x=cbind(mp,summary(dataPrcomp)$importance[2,]+summary(dataPrcomp)$importance[2,1]*0.15),
           labels=paste(round(summary(dataPrcomp)$importance[3,]*100,1),"%",sep=""),
           srt=90)
      # PCA plot
      dev.new()
      plot(cbind(pc1,pc2), pch=21, col="black", bg=colorKey, cex=pc3Size,
           main="PCA", xlab="PC1", ylab="PC2")
      mtext("PC3 (dot size)", side=4)
      # PCA annotaion plot
      dev.new()
      plot(cbind(pc1,pc2), pch=21, cex=pc3Size, col="gray80", bg="gray80", main="PCA",
           xlab="PC1", ylab="PC2", xlim=c(min(pc1)*1.3,max(pc1)*1.3),
           ylim=c(min(pc2)*1.3,max(pc2)*1.3))
      text(cbind(pc1,pc2), labels=names(pc1), cex=0.8)
      # PCA legend
      dev.new()
      plot.new()
      title(main="Legends PCA")
      legend(x=0.05, y=0.9, legend=mainFactors,fill=colorKeyFactors,ncol=length(mainFactors))
      legend(x=0.05, y=0.75, legend=mainFactors,fill=colorKeyFactors)
    }
    
    # Save file:
    if(save == TRUE) {
      dirStat <- dir.create(savedir, recursive=TRUE, showWarnings=FALSE)
      if(dirStat == TRUE) {
        .verb(paste("Creating new directory:",savedir), verbose)
      }
      
      .verb("Saving pca variance plot...", verbose)
      pdfFilePath = paste(savedir,"/","pcaVariance.pdf",sep="")
      if(file.exists(pdfFilePath)) {
        .verb("Warning: pcaVariance.pdf already exists in directory: overwriting old file...", verbose)
      }
      pdf(file=pdfFilePath,paper="a4")
      mp <- barplot(summary(dataPrcomp)$importance[2,],names.arg=colnames(dataPrcomp$rotation),
                    las=3,ylab="Proportion of variance",
                    ylim=c(0,1.25*summary(dataPrcomp)$importance[2,][1]),
                    main="PC imortance")
      text(x=cbind(mp,summary(dataPrcomp)$importance[2,]+summary(dataPrcomp)$importance[2,1]*0.15),
           labels=paste(round(summary(dataPrcomp)$importance[3,]*100,1),"%",sep=""),
           srt=90)
      tmp <- dev.off()
      .verb("...done", verbose)

      .verb("Saving pca plot...", verbose)
      pdfFilePath = paste(savedir,"/","pca.pdf",sep="")
      if(file.exists(pdfFilePath)) {
        .verb("Warning: pca.pdf already exists in directory: overwriting old file...", verbose)
      }
      pdf(file=pdfFilePath,paper="a4")
      plot(cbind(pc1,pc2), pch=21, col="black", bg=colorKey, cex=pc3Size, main="PCA", xlab="PC1", ylab="PC2")
      mtext("PC3 (dot size)", side=4)
      tmp <- dev.off()
      .verb("...done", verbose)
      
      .verb("Saving pca annotation plot...", verbose)
      pdfFilePath = paste(savedir,"/","pcaAnnotation.pdf",sep="")
      if(file.exists(pdfFilePath)) {
        .verb("Warning: pcaAnnotation.pdf already exists in directory: overwriting old file...", verbose)
      }
      pdf(file=pdfFilePath,paper="a4")
      plot(cbind(pc1,pc2), pch=21, cex=pc3Size, col="gray80", bg="gray80", main="PCA",
           xlab="PC1", ylab="PC2", xlim=c(min(pc1)*1.3,max(pc1)*1.3),
           ylim=c(min(pc2)*1.3,max(pc2)*1.3))
      text(cbind(pc1,pc2), labels=names(pc1), cex=0.8)
      tmp <- dev.off()
      .verb("...done", verbose)
      
      .verb("Saving pca legend plot...", verbose)
      pdfFilePath = paste(savedir,"/","pcaLegend.pdf",sep="")
      if(file.exists(pdfFilePath)) {
        .verb("Warning: pcaLegend.pdf already exists in directory: overwriting old file...", verbose)
      }
      pdf(file=pdfFilePath,paper="a4")
      plot.new()
      title(main="Legends")
      legend(x=0.05, y=0.9, legend=mainFactors,fill=colorKeyFactors,ncol=length(mainFactors))
      legend(x=0.05, y=0.75, legend=mainFactors,fill=colorKeyFactors)
      tmp <- dev.off()
      .verb("...done", verbose)
    }
  }
  
  
  # Verbose function:
  .verb <- function(message, verbose) {
    if(verbose == TRUE) {
      message(message)
    }
  }


  # Below is the code that runs the selected functions:
  if(class(arrayData) != "ArrayData") {
    stop("argument arrayData is not of class ArrayData")
  }

  # Run the selected QCs:
    if(rnaDeg == TRUE) {
      if("dataRaw" %in% attributes(arrayData)$names) {
        .runRnaDeg(arrayData, savedir=savedir, save=save)
      } else {
        warning("can not run rnaDeg: argument arrayData does not contain dataRaw")
      }
    }
    if(nuseRle == TRUE) {
      if("dataRaw" %in% attributes(arrayData)$names) {
        .runNuseRle(arrayData, savedir=savedir, save=save)
      } else {
        warning("can not run nuseRle: argument arrayData does not contain dataRaw")
      }
    }
    if(hist == TRUE) {
     .runHist(arrayData, savedir=savedir, save=save)
    }
    if(boxplot == TRUE) {
      .runBoxplot(arrayData, savedir=savedir, save=save)
    }
    if(pca == TRUE) {
      .runPca(arrayData, savedir=savedir, save=save,
              colorFactor=colorFactor, colors=colors)
    }


}