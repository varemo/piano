print.GSAres <- function(x, ...){
   
   cat(paste("Final gene/gene-set association:",x$info$nGenesGSC,"genes and",x$info$nGeneSets,"gene-sets\n"))
   cat("  Details:\n")
   cat(paste("  Calculating gene set statistics from ",x$info$nGenesGSC," out of ",x$info$nGenesStatistics," gene-level statistics\n",sep=""))
   if(x[["signifMethod"]] %in% c("geneSampling","samplePermutation")) cat(paste("  Using all ",x$info$nGenesStatistics," gene-level statistics for significance estimation\n",sep=""))
   cat(paste("  Removed",x$info$removedGenesGSC,"genes from GSC due to lack of matching gene statistics\n"))
   cat(paste("  Removed",x$info$removedGSnoGenes,"gene sets containing no genes after gene removal\n"))
   cat(paste("  Removed additionally",x$info$removedGSsizeLimit,"gene sets not matching the size limits\n"))
   cat(paste("  Loaded additional information for",x$info$nGeneSetsWithAddInfo,"gene sets\n\n"))
   
   tmp <- x[["geneStatType"]]
   if(tmp %in% c("p","p-signed")) tmp <- "p-like"
   else if(tmp %in% c("F","F-signed")) tmp <- "F-like"
   else tmp <- "t-like"
   cat(paste("Gene statistic type: ",tmp,"\n",sep=""))
   cat(paste("Method: ",x[["geneSetStat"]],"\n",sep=""))
   cat(paste("Gene-set statistic name:",x[["gsStatName"]],"\n"))
   if(x[["signifMethod"]] == "geneSampling") tmp <- "Gene sampling"
   if(x[["signifMethod"]] == "samplePermutation") tmp <- "Sample permutation"
   if(x[["signifMethod"]] == "nullDist") tmp <- "Theoretical null distribution"
   cat(paste("Significance: ",tmp,"\n",sep=""))
   cat(paste("Adjustment: ",x[["adjMethod"]],"\n",sep=""))
   cat(paste("Gene set size limit: (",x[["gsSizeLim"]][1],",",x[["gsSizeLim"]][2],")\n",sep=""))
   if(x[["signifMethod"]] != "nullDist" | x[["geneSetStat"]] == "reporter") cat(paste("Permutations:",x[["nPerm"]],"\n"))
   if(x[["geneSetStat"]] == "gsea") cat(paste("GSEA parameter:",x[["gseaParam"]],"\n"))
   cat(paste("Total run time:",round(round(x$runtime[3])/60,2),"min\n"))
   
   
   
   
   
}