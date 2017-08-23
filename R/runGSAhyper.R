#' Gene set analysis with Fisher's exact test
#' 
#' Performs gene set analysis (GSA) based on a list of significant genes and a
#' gene set collection, using Fisher's exact test, returning the gene set
#' p-values.
#' 
#' The statistical test performed is a one-tailed Fisher's exact test on the
#' contingency table with columns "In gene set" and "Not in gene set" and rows
#' "Significant" and "Non-significant" (this is equivalent to a hypergeometric
#' test).
#' 
#' Command run for gene set i:
#' 
#' \code{fisher.test(res$contingencyTable[[i]], alternative="greater")},
#' 
#' the \code{res$contingencyTable} object is available from the object returned
#' from \code{runGSAhyper}.
#' 
#' The main difference between \code{\link{runGSA}} and \code{runGSAhyper} is
#' that \code{runGSA} uses the gene-level statistics (numerical values for each
#' gene) to calculate the gene set p-values, whereas \code{runGSAhyper} only
#' uses the group membership of each gene (in/not in gene set,
#' significant/non-significant). This means that for \code{runGSAhyper} a
#' p-value cut-off for determining significant genes has to be chosen by the
#' user and after this, all significant genes will be seen as equally
#' significant (i.e. the actual p-values are not used). The advantage with
#' \code{runGSAhyper} is that you can use it to find enriched gene sets when
#' you only have a list of interesting genes, without any statistics.
#' 
#' @param genes a vector of all genes in your experiment, or a small list of
#' significant genes.
#' @param pvalues a vector (or object to be coerced into one) of pvalues for
#' genes or a binary vector with 0 for significant genes. Defaults to
#' rep(0,length(genes)), i.e. genes is a vector of genes of interest.
#' @param pcutoff p-value cutoff for significant genes. Defaults to 0 if
#' pvalues are binary. If p-values are spread in [0,1] defaults to 0.05.
#' @param universe a vector of genes that represent the universe. Defaults to
#' genes if pvalues are not all 0. If pvalues are all 0, defaults to all unique
#' genes in gsc.
#' @param gsc a gene set collection given as an object of class \code{GSC} as
#' returned by the \code{\link{loadGSC}} function.
#' @param gsSizeLim a vector of length two, giving the minimum and maximum gene
#' set size (number of member genes) to be kept for the analysis. Defaults to
#' \code{c(1,Inf)}.
#' @param adjMethod the method for adjusting for multiple testing. Can be any
#' of the methods supported by \code{p.adjust}, i.e. \code{"holm"},
#' \code{"hochberg"}, \code{"hommel"}, \code{"bonferroni"}, \code{"BH"},
#' \code{"BY"}, \code{"fdr"} or \code{"none"}.
#' @return A list-like object containing the following elements:
#' 
#' \item{pvalues}{a vector of gene set p-values} \item{p.adj}{a vector of gene
#' set p-values, adjusted for multiple testing} \item{resTab}{a full result
#' table} \item{contingencyTable}{a list of the contingency tables used for
#' each gene set} \item{gsc}{the input gene set collection}
#' @author Leif Varemo \email{piano.rpkg@@gmail.com} and Intawat Nookaew
#' \email{piano.rpkg@@gmail.com}
#' @seealso \pkg{\link{piano}}, \code{\link{loadGSC}}, \code{\link{runGSA}},
#' \code{\link{fisher.test}}, \code{\link{phyper}}, \code{\link{networkPlot}}
#' @examples
#' 
#'    # Load example input data (dummy p-values and gene set collection):
#'    data("gsa_input")
#'    
#'    # Load gene set collection:
#'    gsc <- loadGSC(gsa_input$gsc)
#'    
#'    # Randomly select 100 genes of interest (as an example):
#'    genes <- sample(unique(gsa_input$gsc[,1]),100)
#'       
#'    # Run gene set analysis using Fisher's exact test:
#'    res <- runGSAhyper(genes, gsc=gsc)
#'    
#'    # If you have p-values for the genes and want to make a cutoff for significance:
#'    genes <- names(gsa_input$pvals) # All gene names
#'    p <- gsa_input$pvals # p-values for all genes
#'    res <- runGSAhyper(genes, p, pcutoff=0.001, gsc=gsc)
#'    
#'    # If the 20 first genes are the interesting/significant ones they can be selected
#'    # with a binary vector:
#'    significant <- c(rep(0,20),rep(1,length(genes)-20))
#'    res <- runGSAhyper(genes, significant, gsc=gsc)
#'    
#'    
#' 
runGSAhyper <- function(genes, pvalues, pcutoff, universe, gsc, gsSizeLim=c(1,Inf), adjMethod="fdr") {
   
   ## Argument checking:
   
   # Check gsaSizeLim:
   if(length(gsSizeLim)!=2) stop("argument gsSizeLim should be a vector of length 2")
   
   # Check genes is character vector of unique genes
   if(missing(genes)) {
      stop("argument genes is required")
   } else {
      genes <- as.vector(as.matrix(genes))
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
   
   # Check that there are no genes in GSC not in universe:
   if(!all(unique(unlist(gsc$gsc)) %in% universe)) warning("there are genes in gsc that are not in the universe, these will be removed before analysis")
   
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
   
   # Update gsc, removing genes not in universe, and gene sets not matching size limits:
   gsc <- gsc$gsc
   delInd <- vector()
   for(i in 1:length(gsc)) {
      gs <- gsc[[i]] 
      gs <- gs[gs%in%universe] # In universe!
      if(length(gs) < gsSizeLim[1] | length(gs) > gsSizeLim[2]) delInd <- c(delInd,i)
      gsc[[i]] <- gs
   }
   gsc <- gsc[!c(1:length(gsc))%in%delInd]
   
   # Print info
   message(paste("Analyzing the overrepresentation of ",length(goi)," genes of interest in ",
                 length(gsc)," gene sets, using a background of ",length(bg)," non-interesting genes.",sep=""))
   
   # Loop over genes sets
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
      gs <- gsc[[i]] # In gene set (and in universe!)
      #gs <- gs[gs%in%universe] # In gene set (and in universe!)
      #if(length(gs)!=length(gsc[[i]])) message(paste("Removing ",length(gsc[[i]])-length(gs)," genes from gene set ",i,sep=""))
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
   
   ### NOTE: If changing output structure, check and update networkPlot accordingly!
   
   res <- list()
   res$pvalues <- p
   res$p.adj <- padj
   res$resTab <- resTab
   res$contingencyTable <- contTabList
   res$gsc <- gsc
   return(res)
   
}
