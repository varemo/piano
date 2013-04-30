loadMAdata <- function(datadir=getwd(), setup="setup.txt", dataNorm,
                             platform="NULL", annotation, normalization="plier",
                             filter=TRUE, verbose=TRUE, ...) {

   if(!try(require(affy))) stop("package affy is missing")
   if(!try(require(plier))) stop("package plier is missing")
  
  # Argument check:
  if(!normalization %in% c("plier","rma")) {
    stop("incorrect value of argument normalization")
  }
  if(!platform %in% c("NULL","yeast2")) {
    stop("incorrect value of argument platform")
  }

  # Verbose function:
  .verb <- function(mes, verbose) {
    if(verbose == TRUE) {
      message(mes)
    }
  }
  
  

  # Load the data:
  nCelFiles <- length(dir(path=datadir, pattern = ".*cel", all.files=FALSE,
                          full.names=FALSE, ignore.case = TRUE, recursive=FALSE))
  if(nCelFiles > 0 & missing(dataNorm)) {
    # Load CEL-files
    .verb("Loading CEL files...", verbose)
    dataRaw <- ReadAffy(celfile.path=datadir, ...)
    colnames(exprs(dataRaw)) <- gsub("\\.CEL","",colnames(exprs(dataRaw)), ignore.case=TRUE)
    colnames(exprs(dataRaw)) <- gsub("\\.gz","",colnames(exprs(dataRaw)), ignore.case=TRUE)
    if(sum(duplicated(colnames(exprs(dataRaw)))) > 0) stop("found samples with identical names")
    .verb("...done", verbose)
  } else if(!missing(dataNorm)) {
    if(class(dataNorm) == "character") {
      # If no CEL-files, or if selected, load txt-file
      dataFilePath <- paste(datadir, "/", dataNorm, sep="")
      if(!file.exists(dataFilePath)) {
        stop("could not find the data file")
      }
      .verb("Loading data in text file...", verbose)
      dataNorm <- as.data.frame(read.delim(dataFilePath, header=TRUE, sep="\t",
                                row.names=1, as.is=TRUE,quote=""),stringsAsFactors=FALSE)
      colnames(dataNorm) = gsub("\\.CEL","",colnames(dataNorm), ignore.case=TRUE)
      colnames(dataNorm) = gsub("\\.gz","",colnames(dataNorm), ignore.case=TRUE)
      if(sum(duplicated(colnames(dataNorm))) > 0) stop("found samples with identical names")
      .verb("...done", verbose)
    #
    } else {
     dataNorm <- as.data.frame(dataNorm,stringsAsFactors=FALSE)
    }
  } else {
    stop("could not find any data files in directory")
  }
  # This (above) creates object 'dataRaw' or 'dataNorm' depending on the input
  # (cel or txt).
  
  
  
  # Load the setup:
  if(class(setup) == "character") {
    setupFilePath <- paste(datadir, "/", setup, sep="")
    if(!file.exists(setupFilePath)) {
      stop("could not find the setup file")
    }
    .verb("Loading setup file...", verbose)
    setup <- as.data.frame(read.delim(setupFilePath, header=TRUE, sep="\t",
                           row.names=1, as.is=TRUE,quote=""),stringsAsFactors=FALSE)
    rownames(setup) = gsub("\\.CEL","",rownames(setup), ignore.case=TRUE)
    rownames(setup) = gsub("\\.gz","",rownames(setup), ignore.case=TRUE)
    .verb("...done", verbose)
  } else {
    setup <- as.data.frame(setup, stringsAsFactors=FALSE)
  }
  
  
  # Normalize the raw data:
  if(exists("dataRaw", inherits=FALSE)) {
    # iterplier qubic spline
    if(normalization == "plier") {
      .verb("Preprocessing using PLIER with cubic spline normalization...", verbose)
       dataNorm <- normalize.AffyBatch.qspline(dataRaw, type="pmonly", verbose=FALSE)
       tmp <- suppressWarnings(tmp <- capture.output(dataNorm <- justPlier(dataNorm,normalize=FALSE, 
                                                                           usemm=FALSE, concpenalty=0.08, 
                                                                           plieriteration=30000)))
	    dataNorm <- as.data.frame(exprs(dataNorm),stringsAsFactors=FALSE)
	    colnames(dataNorm) <- gsub("\\.CEL","",colnames(dataNorm), ignore.case=TRUE)
       colnames(dataNorm) <- gsub("\\.gz","",colnames(dataNorm), ignore.case=TRUE)
	    .verb("...done", verbose)
	  } else if(normalization == "rma") {
      .verb("Preprocessing using RMA with quantile normalization...", verbose)
      dataNorm <- rma(dataRaw,verbose=FALSE)
      dataNorm <- as.data.frame(exprs(dataNorm),stringsAsFactors=FALSE)
	    colnames(dataNorm) <- gsub("\\.CEL","",colnames(dataNorm), ignore.case=TRUE)
       colnames(dataNorm) <- gsub("\\.gz","",colnames(dataNorm), ignore.case=TRUE)
	    .verb("...done", verbose)
	  } else if(normalization == "mas5") {
	    .verb("Preprocessing using MAS 5.0 with quantile normalization...", verbose)
      dataNorm <- mas5(dataRaw,verbose=FALSE)
      dataNorm <- as.data.frame(exprs(dataNorm),stringsAsFactors=FALSE)
	    colnames(dataNorm) <- gsub("\\.CEL","",colnames(dataNorm), ignore.case=TRUE)
	    colnames(dataNorm) <- gsub("\\.gz","",colnames(dataNorm), ignore.case=TRUE)
	    .verb("...done", verbose)
	  }
  } else {
    .verb("Text file data: No normalization performed.", verbose)
  }


  # Check annotation
  if(exists("dataRaw", inherits=FALSE)) {
    if(cdfName(dataRaw) == "Yeast_2") {
      platform <- "yeast2"
    }
  }
  if(platform == "yeast2" & missing(annotation)) {
    annotationInfo <- "yeast2"
  } else if(!missing(annotation)) {
    annotationInfo <- "asArgument"
  } else {
    annotationInfo <- "none"
  }
  
  
  if(annotationInfo != "none") {
    if(annotationInfo == "yeast2") {
      if(!try(require(yeast2.db))) stop("package yeast2.db is needed for annotationInfo='yeast2'")
      # Annotate the probes using the yeast2.db package:
      .verb("Creating anotation...", verbose)
      # Gene name
      geneName <- yeast2ORF
    	geneName <- toTable(geneName)
      # Chromosome location
      chromosome <- yeast2CHRLOC
      chromosome <- toTable(chromosome)
      chromosome <- chromosome[,c(1,3,2)]
      # Probe id:s (corresponding to those in dataNorm)
      probeID <- as.data.frame(rownames(dataNorm),stringsAsFactors=FALSE)
    	colnames(probeID) <- "probeID"
    	# Annotation data frame
    	annot <- merge(probeID,geneName,by.x="probeID",by.y="probe_id",all.x=TRUE)
      annot <- merge(annot,chromosome,by.x="probeID",by.y="probe_id",all.x=TRUE)
      rownames(annot) <- annot$probeID
      annot <- annot[2:ncol(annot)]
      colnames(annot) <- c("geneName","chromosome","start") # <- remove sys.name?
      .verb("...done", verbose)
      
    } else if(annotationInfo == "asArgument") {
      # Else annotate from annotation-argument:
      if(class(annotation) == "character") {
        .verb("Creating anotation...", verbose)
        annotFilePath <- paste(datadir, "/", annotation, sep="")
        if(!file.exists(annotFilePath)) {
          stop("could not find the annotation file")
        }
        annot <- as.data.frame(read.delim(annotFilePath, header=TRUE, sep="\t",
                               row.names=1, as.is=TRUE,quote=""),stringsAsFactors=FALSE)
        if(ncol(annot) != 3) {
          stop("provided annotation has to contain 3 columns")
        }
        colnames(annot) <- c("geneName","chromosome","start")
        .verb("...done", verbose)
      } else {
        annot <- as.data.frame(annotation,stringsAsFactors=FALSE)
      }
    }
  } else {
    warning("no annotation created, may cause limitation in downstream functions")
  }
  
  
  if(filter == TRUE & exists("annot", inherits=FALSE)) {
    # Remove unmapped probes:
      .verb("Removing unmapped probes...", verbose)
    	mappedProbes <- rownames(annot)[!is.na(annot[,1])]
    	probes <- rownames(dataNorm)
    	dataNorm <- dataNorm[probes %in% mappedProbes,]
    	annot <- annot[rownames(annot) %in% mappedProbes,]
    	.verb("..done", verbose)
	} else if(filter == TRUE & !exists("annot", inherits=FALSE)) {
	  warning("annotation required for filtering, filtering step is omitted")
	}

  
  
  # Check sample name consistency:
  tmp1 <- length(rownames(setup))
  tmp2 <- length(colnames(dataNorm))
  tmp3 <- length(c(1:tmp1)[rownames(setup) %in% colnames(dataNorm)])
  tmp4 <- length(c(1:tmp2)[colnames(dataNorm) %in% rownames(setup)])
  if(tmp1 != tmp2 | tmp3 != tmp4) {
    stop("inconsistant sample names in dataNorm and setup")
  }

  
  # Construct ArrayData object as return:
  if(exists("dataRaw", inherits=FALSE) & exists("annot", inherits=FALSE)) {
    arrayData <- list(dataRaw=dataRaw, dataNorm=dataNorm, setup=setup, annotation=annot)
  } else if(exists("annot", inherits=FALSE)){
    arrayData <- list(dataNorm=dataNorm, setup=setup, annotation=annot)
  } else if(exists("dataRaw", inherits=FALSE)) {
    arrayData <- list(dataRaw=dataRaw, dataNorm=dataNorm, setup=setup)
  } else {
    arrayData <- list(dataNorm=dataNorm, setup=setup)
  }
  class(arrayData) = "ArrayData"

  return(arrayData)
  
}
