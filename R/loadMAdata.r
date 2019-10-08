#' Load and preprocess microarray data
#'
#' Loads, preprocesses and annotates microarray data to be further used by
#' downstream functions in the \pkg{\link{piano}} package.
#'
#' This function requires at least two inputs: (1) data, either CEL files in
#' the directory specified by \code{datadir} or normalized data specified by
#' \code{dataNorm}, and (2) experimental setup specified by \code{setup}.
#'
#' The setup shold be either a tab delimited text file with column headers or a
#' \code{data.frame}. The first column should contain the names of the CEL
#' files or the column names used for the normalized data, please be sure to
#' use names valid as column names, e.g. avoid names starting with numbers.
#' Additional columns should assign attributes in some category to each array.
#' (For an example run the example below and look at the object
#' \code{myArrayData$setup}.)
#'
#' The \pkg{piano} package is customized for yeast 2.0 arrays and annotation
#' will work automatically, if the cdfName of the arrays equals \emph{Yeast_2}.
#' If using normalized yeast 2.0 data as input, the user needs to set the
#' argument \code{platform="yeast2"} to tell the function to use yeast
#' annotation. If other platforms than yeast 2.0 is used, set
#' \code{platform=NULL} (default) and supply appropriate annotation by the
#' argument \code{annotation}. Note that the cdfName will override
#' \code{platform}, so it can still be set to \code{NULL} for yeast 2.0 CEL
#' files. Note also that \code{annotation} overrides \code{platform}, so if the
#' user wants to use an alternative annotation for yeast, this can be done
#' simply by specifying this in \code{annotation}.
#'
#' The annotation should have the column headers \emph{Gene name},
#' \emph{Chromosome} and \emph{Chromosome location}. The \emph{Gene name} is
#' used in the heatmap in \code{diffExp} and the \emph{Chromosome} and
#' \emph{Chromosome location} is used by the \code{polarPlot}. The rownames (or
#' first column if using a text file) should contain the \emph{probe IDs}. If
#' using a text file the first column should have the header \emph{probeID} or
#' similar. The filtering step discards all probes not listed in the
#' annotation.
#'
#' Normalization is performed on all CEL file data using one of the Affymetrix
#' methods: PLIER (\code{"plier"}) as implemented by
#' \code{\link[plier:justPlier]{justPlier}}, RMA (Robust Multi-Array Average)
#' (\code{"rma"}) expression measure as implemented by
#' \code{\link[affy:rma]{rma}} or MAS 5.0 expression measure \code{"mas5"} as
#' implemented by \code{\link[affy:mas5]{mas5}}.
#'
#' It is possible to pass additional arguments to
#' \code{\link[affy:read.affybatch]{ReadAffy}}, e.g.  \code{cdfname} as this
#' might be required for some types of CEL files.
#'
#' @param datadir character string giving the directory in which to look for
#' the data. Defaults to \code{getwd()}.
#' @param setup character string giving the name of the file containing the
#' experimental setup, or an object of class \code{data.frame} or similar
#' containing the experimental setup. Defaults to \code{"setup.txt"}, see
#' details below for more information.
#' @param dataNorm character string giving the name of the normalized data, or
#' an object of class \code{data.frame} or similar containing the normalized
#' data. Only to be used if the user wishes to start with normalized data
#' rather then CEL files.
#' @param platform character string giving the name of the platform, can be
#' either \code{"yeast2"} or \code{NULL}. See details below for more
#' information.
#' @param annotation character string giving the name of the annotation file,
#' or an object of class \code{data.frame} or similar containing the annotation
#' information. The annotation should consist of the columns \emph{Gene name},
#' \emph{Chromosome} and \emph{Chromosome location}. Not required if
#' \code{platform="yeast2"}.
#' @param normalization character string giving the normalization method, can
#' be either \code{"plier"}, \code{"rma"} or \code{"mas5"}. Defaults to
#' \code{"plier"}.
#' @param filter should the data be filtered? If \code{TRUE} then probes not
#' present in the annotation will be discarded. Defaults to \code{TRUE}.
#' @param verbose verbose? Defaults to \code{TRUE}.
#' @param \dots additional arguments to be passed to \code{ReadAffy}.
#' @return An \code{ArrayData} object (which is essentially a \code{list}) with
#' the following elements:
#'
#' \item{dataRaw}{raw data as an AffyBatch object}
#' \item{dataNorm}{\code{data.frame} containing normalized expression values}
#' \item{setup}{\code{data.frame} containing experimental setup}
#' \item{annotation}{\code{data.frame} containing annotation}
#'
#' Depending on input arguments the \code{ArrayData} object may not include
#' \code{dataRaw} and/or \code{annotation}.
#' @author Leif Varemo \email{piano.rpkg@@gmail.com} and Intawat Nookaew
#' \email{piano.rpkg@@gmail.com}
#' @seealso \pkg{\link{piano}}, \code{\link{runQC}}, \code{\link{diffExp}},
#' \code{\link[affy:read.affybatch]{ReadAffy}},
#' \code{\link[affy:expresso]{expresso}},
#' \code{\link[plier:justPlier]{justPlier}}, \code{\link[yeast2.db:yeast2BASE]{yeast2.db}}
#' @references Gautier, L., Cope, L., Bolstad, B. M., and Irizarry, R. A.  affy
#' - analysis of Affymetrix GeneChip data at the probe level.
#' \emph{Bioinformatics.} \bold{20}, 3, 307-315 (2004).
#' @examples
#'
#'   # Get path to example data and setup files:
#'   dataPath <- system.file("extdata", package="piano")
#'
#'   # Load normalized data:
#'   myArrayData <- loadMAdata(datadir=dataPath, dataNorm="norm_data.txt.gz", platform="yeast2")
#'
#'   # Print to look at details:
#'   myArrayData
#'
#'
loadMAdata <- function(datadir=getwd(), setup="setup.txt", dataNorm,
                             platform="NULL", annotation, normalization="plier",
                             filter=TRUE, verbose=TRUE, ...) {

   #if(!try(require(affy))) stop("package affy is missing") # old, line below is preferred:
   if (!requireNamespace("affy", quietly = TRUE)) stop("package affy is missing")
   #if(!try(require(plier))) stop("package plier is missing") # old, line below is preferred:
   if (!requireNamespace("plier", quietly = TRUE)) stop("package plier is missing")

  # Argument check:
  if(!normalization %in% c("plier","rma","mas5")) {
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
    dataRaw <- affy::ReadAffy(celfile.path=datadir, ...)
    colnames(exprs(dataRaw)) <- gsub("\\.CEL","",colnames(exprs(dataRaw)), ignore.case=TRUE)
    colnames(exprs(dataRaw)) <- gsub("\\.gz","",colnames(exprs(dataRaw)), ignore.case=TRUE)
    if(sum(duplicated(colnames(exprs(dataRaw)))) > 0) stop("found samples with identical names")
    .verb("...done", verbose)
  } else if(!missing(dataNorm)) {
    if(is(dataNorm, "character")) {
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
  if(is(setup, "character")) {
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
       dataNorm <- affy::normalize.AffyBatch.qspline(dataRaw, type="pmonly", verbose=FALSE)
       tmp <- suppressWarnings(tmp <- capture.output(dataNorm <- plier::justPlier(dataNorm,normalize=FALSE,
                                                                           usemm=FALSE, concpenalty=0.08,
                                                                           plieriteration=30000)))
	    dataNorm <- as.data.frame(exprs(dataNorm),stringsAsFactors=FALSE)
	    colnames(dataNorm) <- gsub("\\.CEL","",colnames(dataNorm), ignore.case=TRUE)
       colnames(dataNorm) <- gsub("\\.gz","",colnames(dataNorm), ignore.case=TRUE)
	    .verb("...done", verbose)
	  } else if(normalization == "rma") {
      .verb("Preprocessing using RMA with quantile normalization...", verbose)
      dataNorm <- affy::rma(dataRaw,verbose=FALSE)
      dataNorm <- as.data.frame(exprs(dataNorm),stringsAsFactors=FALSE)
	    colnames(dataNorm) <- gsub("\\.CEL","",colnames(dataNorm), ignore.case=TRUE)
       colnames(dataNorm) <- gsub("\\.gz","",colnames(dataNorm), ignore.case=TRUE)
	    .verb("...done", verbose)
	  } else if(normalization == "mas5") {
	    .verb("Preprocessing using MAS 5.0 with quantile normalization...", verbose)
      dataNorm <- affy::mas5(dataRaw,verbose=FALSE)
      dataNorm <- as.data.frame(log2(exprs(dataNorm)),stringsAsFactors=FALSE)
	    colnames(dataNorm) <- gsub("\\.CEL","",colnames(dataNorm), ignore.case=TRUE)
	    colnames(dataNorm) <- gsub("\\.gz","",colnames(dataNorm), ignore.case=TRUE)
	    .verb("...done", verbose)
	  }
  } else {
    .verb("Text file data: No normalization performed.", verbose)
  }


  # Check annotation
  if(exists("dataRaw", inherits=FALSE)) {
    if(affy::cdfName(dataRaw) == "Yeast_2") {
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
      #if(!try(require(yeast2.db))) stop("package yeast2.db is needed for annotationInfo='yeast2'") # old, line below is preferred:
      if (!requireNamespace("yeast2.db", quietly = TRUE)) stop("package yeast2.db is missing")
      if (!requireNamespace("AnnotationDbi", quietly = TRUE)) stop("package AnnotationDbi is missing")
      # Annotate the probes using the yeast2.db package:
      .verb("Creating annotation...", verbose)
      # Gene name
      geneName <- yeast2.db::yeast2ORF
    	geneName <- AnnotationDbi::toTable(geneName)
      # Chromosome location
      chromosome <- yeast2.db::yeast2CHRLOC
      chromosome <- AnnotationDbi::toTable(chromosome)
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
      if(is(annotation, "character")) {
        .verb("Creating annotation...", verbose)
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
        annot <- annotation
        annot[,1] <- as.character(annotation[,1])
        annot[,2] <- as.character(annotation[,2])
      }

      # Check for NAs:
      suppressWarnings(tmp <- as.numeric(as.character(annot[,3])))
      if(!all(!is.na(tmp))) stop("the chromosome location in annotation has to be numerical")
      annot[,3] <- tmp
      # Remove mappings not in data:
      annot <- annot[rownames(annot)%in%rownames(dataNorm),]
      # Check for duplicates:
      if(length(rownames(annot))!=length(unique(rownames(annot)))) {
         stop("the annotation contains Gene name duplicates")
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
      # Remove mappings not in data:
    	annot <- annot[rownames(annot) %in% rownames(dataNorm),]
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
  class(arrayData) <- "ArrayData"

  return(arrayData)

}
