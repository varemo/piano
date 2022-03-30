#' Perform differential expression analysis
#' 
#' Identifies differentially expressed genes by using the linear model approach
#' of \pkg{limma}. Optionally produces a Venn diagram, heatmap, Polar plot and
#' volcano plot.
#' 
#' This function uses \pkg{\link[limma:01Introduction]{limma}} to calculate p-values measuring
#' differential expression in the given \code{contrasts}. The
#' \code{uniqueFactors} given by \code{\link{extractFactors}} can be used to
#' define a contrast vector, where each element should be a character string on
#' the form \code{"uniqueFactorA - uniqueFactorB"}, note the space surrounding
#' the \code{-}.  (See the example below and for \code{\link{extractFactors}}.)
#' 
#' If appropriate annotation is missing for the \code{ArrayData} object the
#' user can suppply this as \code{chromosomeMapping}. This should be either a
#' \code{data.frame} or a tab delimited text file and include the columns
#' \emph{chromosome} with the chromosome name and \emph{chromosome location}
#' containing the starting position of each gene. A \code{-} sign can be used
#' to denote the antisense strand but this will be disregarded while plotting.
#' The rownames should be \emph{probe IDs} or, if using a text file, the first
#' column with a column header should contain the \emph{probe IDs}.
#' 
#' Note that the \code{fitMethod="robust"} may need longer time to run.
#' 
#' A Venn diagram can be drawn for up to five contrasts (\code{diffExp()} will
#' use \code{vennDiagram}).
#' 
#' The heatmap shows normalized expression values of the genes that pass the
#' \code{heatmapCutoff} in at least one contrast.
#' 
#' A volcano plot is produced for each contrast showing magnitude of change
#' versus significance.
#' 
#' The Polar plot sorts the genes according to chromosomal location, for each
#' chromosome starting with unknown positions followed by increasing number in
#' the \emph{chromosome location} column. Genes which do not map to any
#' chromosome are listed as U for unknown. The radial lines in the Polar plot
#' are -log10 scaled p-values, so that a longer line means a smaller p-value.
#' This gives an overview of the magnitude of differential expression for each
#' contrast.
#' 
#' Typical usages are: \preformatted{ # Identify significantly changed genes in
#' 'm1' and 'm2' compared to 'wt': diffExp(arrayData, contrasts=c("m1 - wt",
#' "m2 - wt")) }
#' 
#' @param arrayData an object of class \code{ArrayData}.
#' @param contrasts a character vector giving the contrasts to be tested for
#' differential expression.  Use \code{\link{extractFactors}} to get allowed
#' contrasts.
#' @param chromosomeMapping character string giving the name of the chromosome
#' mapping file, or an object of class \code{data.frame} or similar containing
#' the chromosome mapping. Required for the Polar plot if the \code{ArrayData}
#' object lacks annotation information. See details below.
#' @param fitMethod character string giving the fitting method used by
#' \code{lmFit}. Can be either \code{"ls"} for least squares (default) or
#' \code{"robust"} for robust regression.
#' @param adjustMethod character string giving the method to use for adjustment
#' of multiple testing.  Can be \code{"holm"}, \code{"hochberg"},
#' \code{"hommel"}, \code{"bonferroni"}, \code{"BH"}, \code{"BY"}, \code{"fdr"}
#' (default) or \code{"none"}. See \code{p.adjust} for details.
#' @param significance number giving the significance cutoff level for the Venn
#' diagram and the horizontal line drawn in the volcano plot.  Defaults to
#' \code{0.001}.
#' @param plot should plots be produced? Set either to \code{TRUE} (default) or
#' \code{FALSE} to control all plots, or to a character vector with any
#' combination of \code{"venn"}, \code{"heatmap"}, \code{"polarplot"} and
#' \code{"volcano"}, to control the single plots (e.g.
#' \code{plot=c("venn","polarplot")} or \code{plot="heatmap"}).
#' @param heatmapCutoff number giving the significance cutoff level for the
#' heatmap. Defaults to \code{1e-10}.
#' @param volcanoFC number giving the x-coordinates of the vertical lines drawn
#' in the volcano plot. Defaults to \code{2}.
#' @param colors character vector of colors to be used by the Venn diagram and
#' Polar plot.
#' @param save should the figures and p-values be saved? Defaults to
#' \code{FALSE}.
#' @param verbose verbose? Defaults to \code{TRUE}.
#' @return A \code{list} with elements:
#' 
#' \item{pValues}{\code{data.frame} containing adjusted p-values (according to
#' argument \code{adjustMethod}) for each contrast}
#' \item{foldChanges}{\code{data.frame} containing log2 fold changes for each
#' contrast} \item{resTable}{a \code{list} with an element for each contrast,
#' each being a \code{data.frame} with full result information}
#' \item{vennMembers}{\code{list} containing the gene members of each area of
#' the Venn diagram (only returned when a Venn diagram is drawn)}
#' @author Leif Varemo \email{piano.rpkg@@gmail.com} and Intawat Nookaew
#' \email{piano.rpkg@@gmail.com}
#' @seealso \pkg{\link{piano}}, \code{\link{loadMAdata}},
#' \code{\link{extractFactors}}, \code{\link{polarPlot}}, \code{\link{runGSA}},
#' \pkg{\link[limma:01Introduction]{limma}}, \code{\link{venn}}, \code{\link{heatmap.2}}
#' @references Smyth, G. K. (2005). Limma: linear models for microarray data.
#' In: 'Bioinformatics and Computational Biology Solutions using R and
#' Bioconductor'. R. Gentleman, V. Carey, S. Dudoit, R. Irizarry, W.  Huber
#' (eds), Springer, New York, pages 397--420.
#' @examples
#' 
#'   # Get path to example data and setup files:
#'   dataPath <- system.file("extdata", package="piano")
#' 
#'   # Load normalized data:
#'   myArrayData <- loadMAdata(datadir=dataPath, dataNorm="norm_data.txt.gz", platform="yeast2")
#' 
#'   # Perform differential expression analysis:
#'   pfc <- diffExp(myArrayData, contrasts=c("aerobic_Clim - anaerobic_Clim",
#'                                           "aerobic_Nlim - anaerobic_Nlim"))
#'                                           
#'   # Order the genes according to p-values, for aerobic_Clim vs anaerobic_Clim:
#'   o <- order(pfc$resTable$'aerobic_Clim - anaerobic_Clim'$P.Value)
#'   
#'   # Display statistics for the top 10 significant genes:
#'   pfc$resTable$'aerobic_Clim - anaerobic_Clim'[o[1:10],]
#' 
#' 
diffExp <- function(arrayData, contrasts, chromosomeMapping,
                                   fitMethod="ls", adjustMethod="fdr", significance=0.001, 
                                   plot=TRUE, heatmapCutoff=1e-10, volcanoFC=2,
                                   colors=c("red","green","blue","yellow","orange",
                                   "purple","tan","cyan","gray60","black"),
                                   save=FALSE, verbose=TRUE) {

  #if(!try(require(limma))) stop("package limma is missing") # old, below is preferred:
  if (!requireNamespace("limma", quietly = TRUE)) stop("package limma is missing")
  
  # Argument check:
  if(missing(contrasts)) stop("argument contrasts is not defined")
  if(!fitMethod %in% c("ls","robust")) {
    stop("incorrect value of argument fitMethod")
  }
  if(!adjustMethod %in% c("holm","hochberg","hommel","bonferroni","BH","BY","fdr","none")) {
    stop("incorrect value of argument adjustMethod")
  }
  if(is(plot, "logical")) {
     if(plot) {
        venn <- heatmap <- polarPlot <- volcano <- TRUE  
     } else {
        venn <- heatmap <- polarPlot <- volcano <- FALSE
     }
  } else if(is(plot, "character")) {
     venn <- heatmap <- polarPlot <- volcano <- FALSE
     if("venn" %in% plot) venn <- TRUE
     if("heatmap" %in% plot) heatmap <- TRUE
     if("polarplot" %in% plot) polarPlot <- TRUE
     if("volcano" %in% plot) volcano <- TRUE
  } else {
     stop("argument plot has to be either TRUE, FALSE or a character string")
  }
  
  savedirFig <- paste(getwd(),"/Piano_Results/Figures/DifferentialExpression", sep="")
  savedirPval <- paste(getwd(),"/Piano_Results/pValues/genes",sep="")
  
  # Verbose function:
  .verb <- function(message, verbose) {
    if(verbose == TRUE) {
      message(message)
    }
  }
  
  saveFig <- save
  savePval <- save
  
  
  # Get factors from setup and sort according to dataNorm columns
  factors <- extractFactors(arrayData)


  # Run lmFit
  .verb("Fitting linear models...", verbose)
  factors <- factor(factors$factors[,1])
  designMatrix <- model.matrix(~0+factors)
  colnames(designMatrix) <- levels(factors)
  dataForLimma <- arrayData$dataNorm
  fitLm <- limma::lmFit(dataForLimma, design=designMatrix, method=fitMethod, maxit=200)


  # Run ebayes
  contrastMatrix <- limma::makeContrasts(contrasts=contrasts, levels=levels(factors))
  fitContrasts <- limma::contrasts.fit(fitLm,contrasts=contrastMatrix)
  fitContrasts <- limma::eBayes(fitContrasts)
  .verb("...done", verbose)
  
  
  # Venn diagram
  if(venn == TRUE) {
    if(length(contrasts) <= 5) {
      .verb("Generating Venn diagrams...", verbose)
      vennInfo <- limma::decideTests(fitContrasts,adjust.method=adjustMethod,p.value=significance)
      # Plot
      if(saveFig == FALSE) {
        #if(length(contrasts) <= 3) {
          dev.new()
          limma::vennDiagram(vennInfo, cex=1, main=paste("Venn diagram (p-value adjustment: ",adjustMethod,", p<",significance,")",sep=""),
                      circle.col=colors, names=LETTERS[1:length(contrasts)])
          dev.new()
          plot.new()
          title(main="Legend Venn diagram")
          #legend(x=0.05, y=0.8, legend=colnames(fitContrasts),fill=colors)
          for(i in 1:length(contrasts)) {
             text(x=0.2, y=(1-0.05*i), paste(LETTERS[i],": ",contrasts[i],sep=""), adj=c(0,0))
          }

        #} else {
          #dev.new()
          #vennInfoTmp <- as.data.frame(vennInfo,stringsAsFactors=FALSE)
          #vennInfoTmp[vennInfoTmp == -1] <- 1
          #colnames(vennInfoTmp) <- LETTERS[1:ncol(vennInfoTmp)]
          #venn(vennInfoTmp)
          
          #dev.new()
          #plot.new()
          #title(main="Legend Venn diagram")
          #for(i in 1:length(contrasts)) {
          #  text(x=0.2, y=(1-0.05*i), paste(LETTERS[i],": ",contrasts[i],sep=""), adj=c(0,0))
          #}
        #}
      }
      .verb("...done", verbose)
      #Save
      if(saveFig == TRUE) {
        dirStat <- dir.create(savedirFig, recursive=TRUE, showWarnings=FALSE)
        if(dirStat == TRUE) {
          .verb(paste("Creating new directory:",savedirFig), verbose)
        }
        # Venn diagram
        vennFileName = paste("venn_pval_adj_",adjustMethod,".pdf",sep="")
        vennFilePath = paste(savedirFig,"/",vennFileName,sep="")
        .verb("Saving Venn diagram...", verbose)
        if(file.exists(vennFilePath)) {
          .verb(paste("Warning: ",vennFileName," already exists in directory: overwriting old file...",sep=""), verbose)
        }
        #if(length(contrasts) <= 3) {
        pdf(file=vennFilePath,paper="a4")
        limma::vennDiagram(vennInfo, cex=1, main=paste("Venn diagram (p-value adjustment: ",adjustMethod,", p<",significance,")",sep=""),
                      circle.col=colors, names=LETTERS[1:length(contrasts)])
        #} else {
          #vennInfoTmp <- as.data.frame(vennInfo,stringsAsFactors=FALSE)
          #vennInfoTmp[vennInfoTmp == -1] <- 1
          #colnames(vennInfoTmp) <- LETTERS[1:ncol(vennInfoTmp)]
          #pdf(file=vennFilePath, paper="a4", onefile=FALSE)
          #venn(vennInfoTmp)
        #}
        tmp <- dev.off()
        .verb("...done", verbose)
        
        # Venn legend
        vennFileName = "vennLegend.pdf"
        vennFilePath = paste(savedirFig,"/",vennFileName,sep="")
        .verb("Saving Venn legend...", verbose)
        if(file.exists(vennFilePath)) {
          .verb("Warning: vennLegend.pdf already exists in directory: overwriting old file...", verbose)
        }
        pdf(file=vennFilePath, paper="a4")
        plot.new()
        #if(length(contrasts) <= 3) {
          #title(main="Legend Venn diagram")
          #legend(x=0.05, y=0.8, legend=colnames(fitContrasts),fill=c("red","green","blue"))
        #} else {
          title(main="Legend Venn diagram")
          for(i in 1:length(contrasts)) {
            text(x=0.2, y=(1-0.05*i), paste(LETTERS[i],": ",contrasts[i],sep=""), adj=c(0,0))
          }
        #}
        tmp <- dev.off()
        .verb("...done", verbose)
      }
    } else {
      warning("can not create a Venn diagram for more than 5 samples")
    }
  }
  
  # Calculate p-values for each contrast:
  if(savePval == TRUE) {
    # Save pval-files, one for each contrast
    dirStat <- dir.create(savedirPval, recursive=TRUE, showWarnings=FALSE)
    if(dirStat == TRUE) {
      .verb(paste("Creating new directory:",savedirPval), verbose)
    }
  }
  pValues <- NA
  foldChanges <- NA
  topTabList <- list()
  for(i in 1:length(contrasts)) {
    if(savePval == TRUE) {
      .verb(paste("Saving p-values for contrast ",colnames(fitContrasts)[i],"...",sep=""), verbose)
      pvalFileName <- paste(colnames(fitContrasts)[i],"_pval_adj_",adjustMethod,".txt",sep="")
      pvalFilePath <- paste(savedirPval,"/",pvalFileName,sep="")
      if(file.exists(pvalFilePath)) {
        .verb(paste("Warning: ",pvalFileName," already exists in directory: overwriting old file...",sep=""), verbose)
      }
    } else {
      .verb(paste("Calculating p-values for contrast ",colnames(fitContrasts)[i],"...",sep=""), verbose)
    }
    topTab <- limma::topTable(fitContrasts,coef=i,adjust.method=adjustMethod,sort.by="none",number=nrow(dataForLimma))
    
    # For limma version compatability:
    if(!"ID" %in% colnames(topTab)) {
      topTab <- cbind(rownames(topTab),topTab)
      colnames(topTab)[1] <- "ID"
    }
    
    # If gene names exist, add them:
    tmp <- arrayData$annotation$geneName
    if(!is.null(tmp)) {
       names(tmp) <- rownames(arrayData$annotation)
       tmp <- as.data.frame(tmp,stringsAsFactors=FALSE)
       topTab <- merge(x=topTab,y=tmp,by.x="ID",by.y="row.names",all.y=TRUE,sort=FALSE)
       topTab <- topTab[,c("ID","tmp",colnames(topTab)[!colnames(topTab)%in%c("ID","tmp")])]
       colnames(topTab)[1:2] <- c("ProbesetID","GeneName") 
    } else {
       colnames(topTab)[1] <- c("ProbesetID") 
    }
    
    # Check that topTab rows are consistent with arrayData$dataNorm rows:
    if(nrow(topTab)!=nrow(arrayData$dataNorm)) {
       stop("failed to order the genes correctly, possibly an issue with the supplied annotation")
    }
    if(!all(topTab[,1]==rownames(arrayData$dataNorm))) {
      stop("failed to order the genes correctly, possibly an issue with the supplied annotation")  
    }
    
    topTabList[[i]] <- topTab
    pValues <- cbind(pValues,topTab$adj.P.Val)  # P-values for all contrasts, to be further used below
    foldChanges <- cbind(foldChanges,topTab$logFC)  # FC for all contrasts
    if(savePval == TRUE) {
      write.table(topTab,sep="\t",file=pvalFilePath,
                  row.names=FALSE, quote=FALSE)
    }
    .verb("...done", verbose)
  }
  pValues <- as.data.frame(pValues[,2:ncol(pValues)],stringsAsFactors=FALSE)
  rownames(pValues) <- topTab[,1]
  colnames(pValues) <- colnames(fitContrasts)
  foldChanges <- as.data.frame(foldChanges[,2:ncol(foldChanges)],stringsAsFactors=FALSE)
  rownames(foldChanges) <- topTab[,1]
  colnames(foldChanges) <- colnames(fitContrasts)
  
  
  # Heatmap
  if(heatmap == TRUE) {
    .verb("Generating heatmap...", verbose)
    genesInHeatmap <- vector()
    for(i in 1:ncol(pValues)) {
      genesInHeatmap <- c(genesInHeatmap, rownames(pValues)[pValues[,i] < heatmapCutoff])
    }
    if(length(unique(genesInHeatmap)) < 2) {
      warning("No genes were selected, change the heatmapCutoff. Omitting heatmap.")
    } else {
      genesInHeatmap <- unique(genesInHeatmap)
      heatmapMatrix <- as.matrix(arrayData$dataNorm[rownames(arrayData$dataNorm) %in% genesInHeatmap,])
      if("annotation" %in% attributes(arrayData)$names) {
         geneNamesHeatmap <- arrayData$annotation[rownames(arrayData$annotation) %in% genesInHeatmap,]
         for(j in 1:nrow(heatmapMatrix)) {
            tmp <- geneNamesHeatmap[rownames(geneNamesHeatmap) == rownames(heatmapMatrix)[j],1]
            if(length(tmp) == 1) {
               rownames(heatmapMatrix)[j] <- tmp
            }
         }
      }
      if(saveFig == TRUE) {
        dirStat <- dir.create(savedirFig, recursive=TRUE, showWarnings=FALSE)
        if(dirStat == TRUE) {
          .verb(paste("Creating new directory:",savedirFig), verbose)
        }
        figFileName = paste("heatmap_genes.pdf",sep="")
        figFilePath = paste(savedirFig,"/",figFileName,sep="")
        .verb("Saving heatmap...", verbose)
        if(file.exists(figFilePath)) {
          .verb(paste("Warning: ",figFileName," already exists in directory: overwriting old file...",sep=""), verbose)
        }
        pdf(file=figFilePath,paper="a4", width=21, height=29)
      } else {
        dev.new()
      }
      heatmap.2(heatmapMatrix, trace="none", scale="none", margins=c(10,5),
                lwid=c(2,10), lhei=c(2,20), key=FALSE)
      if(saveFig == TRUE) {
        tmp <- dev.off()
        figFileName = paste("heatmap_genes_colorkey.pdf",sep="")
        figFilePath = paste(savedirFig,"/",figFileName,sep="")
        if(file.exists(figFilePath)) {
          .verb(paste("Warning: ",figFileName," already exists in directory: overwriting old file...",sep=""), verbose)
        }
        pdf(file=figFilePath,paper="a4", width=4, height=3)
      } else {
        dev.new()
      }
      maColorBar(seq(from=min(heatmapMatrix),to=max(heatmapMatrix),length.out=100),
                 k=5,main="Color key")
      title(xlab="expression values")
      if(saveFig == TRUE) {
        tmp <- dev.off()
      }
      .verb("...done", verbose)
    }
  }  


  # Polar plot
  if(polarPlot == TRUE) {
    if("annotation" %in% attributes(arrayData)$names & missing(chromosomeMapping)) {
      polarPlot(pValues, chromosomeMapping=arrayData$annotation[,c(2,3)],
                colors=colors, save=saveFig, verbose=verbose)
    } else if(!missing(chromosomeMapping)) {
      polarPlot(pValues, chromosomeMapping=chromosomeMapping,
                colors=colors, save=saveFig, verbose=verbose)
    } else {
      warning("can not run PolarPlot: no chromosome mapping provided")
    }
  }
  
  
  # Volcano plot
  if(volcano == TRUE) {
     .verb("Generating volcano plot...", verbose)
     
     if(saveFig == TRUE) {
        dirStat <- dir.create(savedirFig, recursive=TRUE, showWarnings=FALSE)
        if(dirStat == TRUE) {
           .verb(paste("Creating new directory:",savedirFig), verbose)
        }
        figFileName = paste("volcano_plot.pdf",sep="")
        figFilePath = paste(savedirFig,"/",figFileName,sep="")
        .verb("Saving volcano plot...", verbose)
        if(file.exists(figFilePath)) {
           .verb(paste("Warning: ",figFileName," already exists in directory: overwriting old file...",sep=""), verbose)
        }
        pdf(file=figFilePath,paper="a4")
     }
     
     for(i in 1:ncol(pValues)) {
        if(saveFig == FALSE) dev.new()
        plot(x=foldChanges[,i],y=-log10(pValues[,i]),pch=16,col="gray",
             main=contrasts[i],xlab="log2 fold change",ylab="-log10 p-value",xlim=c(max(abs(foldChanges[,i])),-max(abs(foldChanges[,i]))))
        tmp <- cbind(foldChanges[,i],pValues[,i])
        tmp <- tmp[abs(tmp[,1])>=abs(volcanoFC) & tmp[,2]<=significance,]
        if(length(tmp) > 2) {
           points(x=tmp[,1],y=-log10(tmp[,2]),pch=16,col="black")
        } else if(length(tmp) == 2) {
           points(x=tmp[1],y=-log10(tmp[2]),pch=16,col="black")
        }
        abline(h=-log10(significance),col="red",lty=2)
        abline(v=volcanoFC,col="blue",lty=2)
        abline(v=-volcanoFC,col="blue",lty=2)
     }
     
     if(saveFig == TRUE) {
        tmp <- dev.off()
     }
     .verb("...done", verbose)
  }
  
  
  
  # Get Venn gene membership info:
  if(venn == TRUE) {
     combMat <- expand.grid(lapply(1:ncol(pValues),function(x) 0:1))
     signGeneList <- apply(pValues,2,function(x) which(x < significance))
     if(!is(signGeneList, "list")) {
         tmp <- list()
         for(i in 1:ncol(signGeneList)) {
            tmp[[i]] <- signGeneList[,i]  
         }
         signGeneList <- tmp
     }
     vennList <- list()
     for(i in 2:nrow(combMat)) { # skip the first row with all zeros
        if(ncol(combMat)==1) {
           geneInd <- signGeneList[[1]]
        } else {
           geneInd <- Reduce(intersect, signGeneList[combMat[i,]==1])
           geneInd <- geneInd[!geneInd %in% unique(unlist(signGeneList[combMat[i,]==0]))]
        }
        if(length(topTabList[[1]]$GeneName) > 0) {
           vennList[[i-1]] <- topTabList[[1]][geneInd,1:2]
        } else {
           vennList[[i-1]] <- topTabList[[1]][geneInd,1]
        }
        rownames(vennList[[i-1]]) <- NULL
        names(vennList)[i-1] <- paste("Uniquely in ",paste(LETTERS[1:ncol(pValues)][combMat[i,]==1],collapse=""),sep="")
     }
  }
  
  
  
  # Output
  names(topTabList) <- contrasts
  if(venn == TRUE) {
     return(list(pValues=pValues,foldChanges=foldChanges,resTable=topTabList,vennMembers=vennList))
  } else {
     return(list(pValues=pValues,foldChanges=foldChanges,resTable=topTabList))
  }
}
