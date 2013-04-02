geneSetSummary <- function(gsaRes, geneSet) {
   
   test <- 1 # Which comparison, currently only 1 allowed!
   
   # Error check:
   obj <- gsaRes
   if(class(obj) != "GSAres") stop("argument GSAres has to be of class 'GSAres'")
   if(test > ncol(obj$nGenesTot)) stop("argument test is to large")
   gsc <- obj$gsc
   gsInd <- which(names(gsc) == geneSet)
   if(length(gsInd) == 0) stop("could not find the specified gene set, check the spelling")
   if(length(gsInd) > 1) stop("argument geneSet matches more than one gene set")
   
   # Get all stat info:
   tab <- as.data.frame(c(obj$nGenesTot[gsInd,test], obj$statDistinctDir[gsInd,test], obj$statDistinctDirUp[gsInd,test],
                          obj$pDistinctDirUp[gsInd,test], obj$pAdjDistinctDirUp[gsInd,test], obj$statDistinctDirDn[gsInd,test],
                          obj$pDistinctDirDn[gsInd,test], obj$pAdjDistinctDirDn[gsInd,test], obj$statNonDirectional[gsInd,test],
                          obj$pNonDirectional[gsInd,test], obj$pAdjNonDirectional[gsInd,test], obj$nGenesUp[gsInd,test],
                          obj$statMixedDirUp[gsInd,test], obj$pMixedDirUp[gsInd,test], obj$pAdjMixedDirUp[gsInd,test],
                          obj$nGenesDn[gsInd,test], obj$statMixedDirDn[gsInd,test], obj$pMixedDirDn[gsInd,test],
                          obj$pAdjMixedDirDn[gsInd,test]), stringsAsFactors=FALSE)
   
   names <- c("Genes (tot)", "Stat (dist.dir.)", "Stat (dist.dir.up)", "p (dist.dir.up)", "p adj (dist.dir.up)", "Stat (dist.dir.dn)",
              "p (dist.dir.dn)", "p adj (dist.dir.dn)","Stat (non-dir)", "p (non-dir)", "p adj (non-dir)", "Genes (up)",
              "Stat (mix.dir.up)", "p (mix.dir.up)", "p adj (mix.dir.up)", "Genes (dn)", "Stat (mix.dir.dn)",
              "p (mix.dir.dn)", "p adj (mix.dir.dn)")
   
   tab[,2] <- names
   tab <- tab[,c(2,1)]
   colnames(tab) <- c("Name","Value")
   tab <- tab[!is.na(tab[,2]),]
   rownames(tab) <- 1:nrow(tab)
   
   # Gene-set name
   name <- names(gsc)[gsInd]
   
   # Gene-set genes:
   genes <- gsc[[gsInd]]
   tmp <- obj$geneLevelStats
   geneLevelStats <- tmp[rownames(tmp)%in%genes,test]
   tmp <- obj$directions
   if(class(tmp) == "character") {
      directions <- tmp
   } else {
      directions <- tmp[rownames(tmp)%in%genes,test]
   }
   
   # Contrast name:
   #contrastName <- obj$contrastName[test]
   
   # Return result:
   res <- list()
   res$name <- name
   res$geneLevelStats <- geneLevelStats
   res$directions <- directions
   #res$contrastName <- contrastName
   res$stats <- tab
   return(res)
   
}