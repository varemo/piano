#' Gene set analysis summary table
#' 
#' Displays or saves a summary table of the results from \code{\link{runGSA}}.
#' 
#' The table is by default saved as an .xls file, if \code{file} is unused.
#' 
#' @param gsaRes an object of class \code{GSAres}, as returned from
#' \code{runGSA()}.
#' @param save a logical, whether or not to save the table.
#' @param file a character string giving the file name to save to.
#' @return The summary table as a data.frame (returned invisibly if
#' \code{save=TRUE}).
#' @author Leif Varemo \email{piano.rpkg@@gmail.com} and Intawat Nookaew
#' \email{piano.rpkg@@gmail.com}
#' @seealso \pkg{\link{piano}}, \code{\link{runGSA}},
#' \code{\link{networkPlot}}, \code{\link{GSAheatmap}}
#' @examples
#' 
#'    # Load example input data to GSA:
#'    data("gsa_input")
#'    
#'    # Load gene set collection:
#'    gsc <- loadGSC(gsa_input$gsc)
#'       
#'    # Run gene set analysis:
#'    gsares <- runGSA(geneLevelStats=gsa_input$pvals , directions=gsa_input$directions, 
#'                     gsc=gsc, nPerm=500)
#'       
#'    # Summary table:
#'    GSAsummaryTable(gsares)  
#' 
GSAsummaryTable <- function(gsaRes, save=FALSE, file=NULL) {
   
   test <- 1 # which comparison, currently only 1 allowed!
   obj <- gsaRes
   if(!is(obj, "GSAres")) stop("argument GSAres has to be of class 'GSAres'")
   if(!is.null(file)) {
      if(!is(file, "character")) stop("argument file has to be of class 'character'")  
   }
   if(test > ncol(obj$nGenesTot)) stop("argument test is to large")
   
   tab <- data.frame(names(obj$gsc), 
                     cbind(obj$nGenesTot)[,test], 
                     cbind(obj$statDistinctDir)[,test], 
                     cbind(obj$statDistinctDirUp)[,test], 
                     cbind(obj$pDistinctDirUp)[,test],
                     cbind(obj$pAdjDistinctDirUp)[,test], 
                     cbind(obj$statDistinctDirDn)[,test], 
                     cbind(obj$pDistinctDirDn)[,test], 
                     cbind(obj$pAdjDistinctDirDn)[,test],
                     cbind(obj$statNonDirectional)[,test], 
                     cbind(obj$pNonDirectional)[,test], 
                     cbind(obj$pAdjNonDirectional)[,test], 
                     cbind(obj$nGenesUp)[,test],
                     cbind(obj$statMixedDirUp)[,test], 
                     cbind(obj$pMixedDirUp)[,test], 
                     cbind(obj$pAdjMixedDirUp)[,test], 
                     cbind(obj$nGenesDn)[,test],
                     cbind(obj$statMixedDirDn)[,test], 
                     cbind(obj$pMixedDirDn)[,test], 
                     cbind(obj$pAdjMixedDirDn)[,test], stringsAsFactors=F)
   
   colnames(tab) <- c("Name", "Genes (tot)", "Stat (dist.dir)", "Stat (dist.dir.up)", "p (dist.dir.up)", "p adj (dist.dir.up)",
                      "Stat (dist.dir.dn)", "p (dist.dir.dn)", "p adj (dist.dir.dn)", "Stat (non-dir.)", "p (non-dir.)", "p adj (non-dir.)",
                      "Genes (up)", "Stat (mix.dir.up)", "p (mix.dir.up)", "p adj (mix.dir.up)", "Genes (down)",
                      "Stat (mix.dir.dn)", "p (mix.dir.dn)", "p adj (mix.dir.dn)")
   rownames(tab) <- NULL

   # Remove unused columns:
   tab <- tab[,apply(tab,2,function(x) sum(is.na(x)) < nrow(tab))]
   
   if(save) {
      ind <- ""
      while(is.null(file)) {
         tmp <- length(grep(paste("^GSAtab",as.character(ind),".xls$",sep=""),dir()))
         if(tmp > 0) {
            if(ind == "") ind <- 0
            ind <- ind + 1   
         } else {
            file <- paste("GSAtab",as.character(ind),".xls",sep="")
         }
         
      }
      write.table(tab,sep="\t", file=file, row.names=FALSE, quote=FALSE)
   }
   
   if(save) {
      invisible(tab)
   } else {
      return(tab)  
   }
   
}
