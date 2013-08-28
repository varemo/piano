runGSAhyper <- function(genes, pvalues, pcutoff, universe, gsc, adjMethod="fdr") {
   
   ## Argument checking:
   
   # Check genes is character vector of unique genes
   if(missing(genes)) {
      stop("argument genes is required")
   } else {
      if(class(genes) != "character") stop("argument genes should be a character vector")  
      if(length(unique(genes)) != length(genes)) stop("argument genes should contain no duplicated entries")
   }
   
   # If pvalues not set, set all 0
   if(missing(pvalues)) {
      pvalues <- rep(0,length(genes))
   } else {
   
   # If pvalues set, check in [0,1] and numeric vector of length length(genes)
      pvalues <- as.vector(as.matrix(pvalues))
      if(class(pvalues) != "numeric") stop("argument pvalues should be a numeric vector")
      if(length(pvalues) != length(genes)) stop("argument pvalues should be the same length as argument genes")
      if(max(pvalues)>1 | min(pvalues)<0) stop("pvalues need to lie between 0 and 1")
   }
      
   # If pcutoff not set & pvalues binary, set to 0
   if(missing(pcutoff)) {
      if(all(pvalues%in%c(0,1))) {
         pcutoff <- 0
      } else {
   
   # If pcutoff not set & pvalues spread in [0,1], set to 0.05
         pcutoff <- 0.05
      }
   } else {
      
   # If pcutoff set, check length 1 and numeric and in [0,1]
      if(length(pcutoff) != 1 & class(pcutoff) != "numeric") stop("argument pcutoff should be a numeric of length 1")
      if(max(pcutoff)>1 | min(pcutoff)<0) stop("argument pcutoff needs to lie between 0 and 1")
   }
   
   # Check gsc
   if(missing(gsc)) {
      stop("argument gsc needs to be given")  
   } else {
      if(class(gsc) != "GSC") stop("argument gsc should be of class GSC, as returned by the loadGSC function")
   }
   
   # If pvalues are not all 0, and universe not set, set to genes
   if(missing(universe)) {
      if(!all(pvalues==0)) {
         universe <- genes
         message("Using all genes in argument genes as universe.")
      } else {
         
   # If pvalues are all 0, and universe not set, set to unique(unlist(gsc$gsc))
         universe <- unique(unlist(gsc$gsc))
         message("Using all genes present in argument gsc as universe.")
      }
   } else {
      
   # If universe set, check character vector
      if(class(universe) != "character") stop("argument universe should be a character vector")
      if(!all(pvalues==0)) stop("if universe is given, genes should be only the genes of interest, i.e. pvalues should all be set to 0.")
   }
   
   # Check that all(genes%in%universe)==T.
   if(!all(genes%in%universe)) {
      warning("not all genes given by argument genes are present in universe, these will be added to universe")  
      universe <- c(universe,genes[!genes%in%universe])
   }
   
   # Check no duplicates in universe
   if(length(unique(universe))!=length(universe)) stop("argument universe should contain no duplicated entries")
   
   # Check adjMethod:
   tmp <- try(adjMethod <- match.arg(adjMethod, c("holm", "hochberg", "hommel", "bonferroni", 
                                                  "BH", "BY","fdr", "none"), several.ok=FALSE), silent=TRUE)
   if(class(tmp) == "try-error") {
      stop("argument adjMethod set to unknown method")
   }
   
   ## Argument checking done!
   
   
   # Set pvalues==0 to -1e10, this so that p<cutoff will be p<=cutoff if cutoff=0.
   pvalues[pvalues==0] <- -1e-10
   
   # Get genes of interest
   goi <- genes[pvalues<pcutoff] # Significant
   if(length(goi)<1) stop("no genes selected due to too strict pcutoff")
   
   # Get background
   bg <- universe[!universe%in%goi] # Not significant
   
   # Print info
   message(paste("Analyzing the overrepresentation of ",length(goi)," genes of interest in ",
                 length(gsc$gsc)," gene sets, using a background of ",length(bg)," non-interesting genes.",sep=""))
   
   # Loop over genes sets
   gsc <- gsc$gsc
   p <- rep(NA,length(gsc))
   names(p) <- names(gsc)
   padj <- rep(NA,length(gsc))
   names(padj) <- names(gsc)
   contTabList <- list() 
   resTab <- matrix(nrow=length(gsc),ncol=6)
   colnames(resTab) <- c("p-value","Adjusted p-value","Significant (in gene set)",
                         "Non-significant (in gene set)","Significant (not in gene set)",
                         "Non-significant (not in gene set)")
   rownames(resTab) <- names(gsc)
   for(i in 1:length(gsc)) {
      gs <- gsc[[i]] # In gene set
      nogs <- universe[!universe%in%gs] # Not in gene set
      
      
      ctab <- rbind(c(sum(goi%in%gs),sum(goi%in%nogs)),
                    c(sum(bg%in%gs),sum(bg%in%nogs)))
      p[i] <- fisher.test(ctab,alternative="greater")$p.value
      
      #p2[i] <- 1-phyper(sum(goi%in%gs)-1,length(goi),length(bg),length(gs))
      
      rownames(ctab) <- c("Significant","Non-significant")
      colnames(ctab) <- c("Genes in gene set","Genes not in gene set")
      contTabList[[i]] <- ctab
      
      
      
      resTab[i,] <- c(p[i],NA,sum(goi%in%gs),sum(bg%in%gs),sum(goi%in%nogs),sum(bg%in%nogs)) 
      
   }
   
   padj <- p.adjust(p, method=adjMethod)
   resTab[,2] <- padj
   
   res <- list()
   res$pvalues <- p
   res$p.adj <- padj
   res$resTab <- resTab
   res$contingencyTable <- contTabList
   #res$gsc <- gsc
   return(res)
   
}