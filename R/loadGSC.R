loadGSC <- function(file, type="auto", addInfo) {

   # Initial argument checks:
   if(missing(addInfo)) {
      addUserInfo <- "skip"
      addInfo <- "none"
   } else {
      addUserInfo <- "yes"
   }
   
   tmp <- try(type <- match.arg(type, c("auto","gmt","sbml","sif","data.frame"), several.ok=FALSE), silent=TRUE)
   if(class(tmp) == "try-error") {
      stop("argument type set to unknown value")
   }
   
   # Check file extension if type="auto":
   if(type == "auto") {
      if(class(file) == "character") {
         tmp <- unlist(strsplit(file,"\\."))
         type <- tolower(tmp[length(tmp)])
         if(!type %in% c("gmt","sif","sbml","xml")) stop(paste("can not handle .",type," file extension, read manually using e.g. read.delim() and load as data.frame",sep=""))
      } else {
         type <- "data.frame"
      }
   }
   
   
   #************************
   # GMT
   #************************
   
   # Read gmt-file:
   if(type == "gmt") {
      
      con <- file(file)
      tmp <- try(suppressWarnings(open(con)), silent=TRUE)
      if(class(tmp) == "try-error") stop("file could not be read")
      if(addUserInfo == "skip") addInfo <- vector()
      gscList <- list()
      i <- 1
      tmp <- try(suppressWarnings(
      while(length(l<-scan(con,nlines=1,what="character",quiet=T)) > 0) {
         if(addUserInfo == "skip") addInfo <- rbind(addInfo,l[1:2])
         tmp <- l[3:length(l)]
         gscList[[l[1]]] <- unique(tmp[tmp != "" & tmp != " " & !is.na(tmp)])
         i <- i + 1
      }
      ), silent=TRUE)
      if(class(tmp) == "try-error") stop("file could not be read")
      close(con)
      
      # Remove duplicate gene sets:
      gsc <- gscList[!duplicated(names(gscList))]
      if(addUserInfo == "skip") addInfo <- unique(addInfo)
      #info$redundantGS <- length(gscList) - length(gsc)
   
      
   #************************
   # SBML
   #************************
      
   } else if(type %in% c("sbml","xml")) {
      
      require(rsbml)
      # Read sbml file:
      tmp <- try(sbml <- rsbml_read(file))
      if(class(tmp) == "try-error") {
         stop("file could not be read by rsbml_read()")
      }
      
      # Create gsc object:
      gsc <- list()
      for(iReaction in 1:length(reactions(model(sbml)))) {
         
         # Species ID for metabolites in current reaction:
         metIDs <- names(c(reactants(reactions(model(sbml))[[iReaction]]),
                                products(reactions(model(sbml))[[iReaction]])))
         
         # Get gene id:s for genes associated with current reaction:
         geneIDs <- names(modifiers(reactions(model(sbml))[[iReaction]]))
         
         # If any genes found:
         if(length(geneIDs) > 0) {
            
            # Get gene names:
            geneNames <- rep(NA,length(geneIDs))
            for(iGene in 1:length(geneIDs)) {
               geneNames[iGene] <- name(species(model(sbml))[[geneIDs[iGene]]])
            }
            
            # Loop over metabolites for current reaction, add gene names:
            for(iMet in 1:length(metIDs)) {
               gsc[[metIDs[iMet]]] <- c(gsc[[metIDs[iMet]]], geneNames)
            }
         }
      }
      
      # Fix the gene-set names to metabolite names (in place of ids):
      if(length(gsc) == 0) {
         stop("no gene association found")
      } else {
         for(iMet in 1:length(gsc)) {         
            tmp1 <- name(species(model(sbml))[[names(gsc)[iMet]]])
            tmp2 <- compartment(species(model(sbml))[[names(gsc)[iMet]]])
            names(gsc)[iMet] <- paste(tmp1," (",tmp2,")",sep="")
         }
      }
      
      
   #************************
   # SIF
   #************************ 
      
   } else if(type == "sif") {
      tmp <- try(gsc <- as.data.frame(read.delim(file, header=FALSE, quote="", as.is=TRUE), 
                                      stringsAsFactors=FALSE), silent=TRUE)
      if(class(tmp) == "try-error") {
         stop("argument file could not be read and converted into a data.frame")
      }
      
      # Check gsc for three columns:
      if(ncol(gsc)!=3) {
         stop("sif file should contain three columns")  
      }
      
      # Get gsc and addInfo part:
      if(addUserInfo == "skip") addInfo <- gsc[,c(1,2)]
      gsc <- gsc[,c(3,1)]
      
      # Remove redundant rows:
      tmp <- nrow(gsc)
      gsc <- unique(gsc)
      #info$redundantGS <- tmp - nrow(gsc)
      
      # Convert to list object:
      geneSets <- unique(gsc[,2])
      gscList <- list()
      for(iGeneSet in 1:length(geneSets)) {
         gscList[[iGeneSet]] <- gsc[gsc[,2] == geneSets[iGeneSet],1]
      }
      names(gscList) <- geneSets
      gsc <- gscList
      
      
   #************************
   # Data.frame
   #************************
   
   # Gene set collection as data.frame:
   } else if(type == "data.frame") {
      tmp <- try(gsc <- as.data.frame(file, stringsAsFactors=FALSE), silent=TRUE)
      if(class(tmp) == "try-error") {
         stop("argument file could not be converted into a data.frame")
      }
      # Get rid of factors:
      for(i in 1:ncol(gsc)) {
         gsc[,i] <- as.character(gsc[,i])
      }
      
      # Check gsc for two columns:
      if(ncol(gsc)!=2) {
         stop("argument file has to contain exactly two columns")  
      }
      
      # Remove redundant rows:
      tmp <- nrow(gsc)
      gsc <- unique(gsc)
      #info$redundantGS <- tmp - nrow(gsc)
      
      # Convert to list object:
      geneSets <- unique(gsc[,2])
      gscList <- list()
      for(iGeneSet in 1:length(geneSets)) {
         gscList[[iGeneSet]] <- gsc[gsc[,2] == geneSets[iGeneSet],1]
      }
      names(gscList) <- geneSets
      gsc <- gscList
   }
   
   
   #***************************
   # AddInfo
   #***************************
      
   # Additional info as data.frame:
   if(addUserInfo == "yes") {
      tmp <- try(addInfo <- as.data.frame(addInfo, stringsAsFactors=FALSE), silent=TRUE)
      if(class(tmp) == "try-error") {
         stop("failed to convert additional info in argument 'addInfo' into a data.frame")
      }
   }
   
   if(class(addInfo) == "data.frame") {
      
      # Check for 2 columns:
      if(ncol(addInfo) != 2) stop("additional info in argument 'file' or 'addInfo' has to contain 2 columns")
      
      # Check addInfo correlation to gsc:
      tmp <- nrow(addInfo)
      addInfo <- unique(addInfo[addInfo[,1] %in% names(gsc),])
      #info$unmatchedAddInfo <- tmp - nrow(addInfo)
   } else {
      #info$unmatchedAddInfo <- 0     
   }
   
   #********************************
   # Return values:
   #********************************
   
   res <- list(gsc,addInfo)
   names(res) <- c("gsc","addInfo")
   class(res) <- "GSC"
   return(res)
   
}