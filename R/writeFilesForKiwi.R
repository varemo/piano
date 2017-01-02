writeFilesForKiwi <- function(gsaRes, label="", overwrite=FALSE) {
   
   # Error checks:
   listExists <- FALSE
   if(class(gsaRes)=="list") {
      gsaResList <- gsaRes
      gsaRes <- gsaRes[[1]]
      listExists <- TRUE
   }
   if(class(gsaRes)!="GSAres") stop("argument gsaRes should be of class GSAres or a list of GSAres objects")
   if(class(label)!="character" | length(label)!=1) stop("argument label should be a character string")
   if(nchar(label)>0) {
      label <- paste("_",label,sep="")
   }
   
   # Gene set - gene file:
   gsc <- gsaRes$gsc
   gscTable <- vector()
   for(i in 1:length(gsc)) {
      gscTable <- rbind(gscTable,cbind(names(gsc)[i],gsc[[i]]))
   }
   filename <- paste("GSC",label,".txt",sep="")
   if(!file.exists(filename)) {
      write.table(gscTable, file=filename, sep="\t", 
                  col.names=FALSE, row.names=FALSE, quote=FALSE)
   } else {
      if(overwrite) {
         write.table(gscTable, file=filename, sep="\t", 
                     col.names=FALSE, row.names=FALSE, quote=FALSE)
         warning(paste("file ",filename)," was overwritten")
      } else {
         stop("files with identical names already exist, replace names or enforce overwriting")
      }
   }
   
   # Gene-level statistics file:
   if(gsaRes$geneStatType == "t") {
      warning("Only gene-level p-values are currently supported to be exported to Kiwi. Skipping GLS-file.") # Not supported in Kiwi yet
      #t <- gsaRes$geneLevelStats
      #fc <- sign(t)
      #glTable <- merge(p,fc,by=0)
      #colnames(glTable) <- c("g","t","FC")
   } else if(gsaRes$geneStatType == "p-signed") {
      p <- gsaRes$geneLevelStats
      fc <- gsaRes$directions
      glTable <- data.frame(g=rownames(p),p=p[,1],FC=fc[,1],stringsAsFactors=F)
      colnames(glTable) <- c("g","p","FC")
   } else {
     warning("Only gene-level p-values are currently supported to be exported to Kiwi. Skipping GLS-file.") # Not supported in Kiwi yet
   }
   filename <- paste("GLS",label,".txt",sep="")
   if(!file.exists(filename)) {
      write.table(glTable, file=filename, sep="\t", 
                  col.names=TRUE, row.names=FALSE, quote=FALSE)
   } else {
      if(overwrite) {
         write.table(glTable, file=filename, sep="\t", 
                     col.names=TRUE, row.names=FALSE, quote=FALSE)
         warning(paste("file ",filename)," was overwritten")
      } else {
         stop("files with identical names already exist, replace names or enforce overwriting")
      }
   }
   
   # Gene set statistics file:
   if(listExists) {
      suppressWarnings(tmp1 <- consensusHeatmap(gsaResList,cutoff=Inf, adjusted=TRUE, plot=FALSE)$pMat)
      suppressWarnings(tmp2 <- consensusHeatmap(gsaResList,cutoff=Inf, adjusted=FALSE, plot=FALSE)$pMat)
      gsTable <- merge(tmp1,tmp2,by=0)
      colnames(gsTable) <- c("Name",
                             "p adj (dist.dir.dn)","p adj (mix.dir.dn)","p adj (non-dir.)",
                             "p adj (mix.dir.up)","p adj (dist.dir.up)",
                             "p (dist.dir.dn)","p (mix.dir.dn)","p (non-dir.)",
                             "p (mix.dir.up)","p (dist.dir.up)")

      filename <- paste("GSS",label,".txt",sep="")
      if(!file.exists(filename)) {
         write.table(gsTable, file=filename, sep="\t", 
                     col.names=TRUE, row.names=FALSE, quote=FALSE)
      } else {
         if(overwrite) {
            write.table(gsTable, file=filename, sep="\t", 
                        col.names=TRUE, row.names=FALSE, quote=FALSE)
            warning(paste("file ",filename)," was overwritten")
         } else {
            stop("files with identical names already exist, replace names or enforce overwriting")
         }
      }
   } else {
      filename <- paste("GSS",label,".txt",sep="")
      if(!file.exists(filename)) {
         GSAsummaryTable(gsaRes,save=TRUE,file=filename)
      } else {
         if(overwrite) {
            GSAsummaryTable(gsaRes,save=TRUE,file=filename)
            warning(paste("file ",filename)," was overwritten")
         } else {
            stop("files with identical names already exist, replace names or enforce overwriting")
         }
      }
   }
}