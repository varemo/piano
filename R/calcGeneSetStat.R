calcGeneSetStat <- function(selectedStats, method, statistics=NULL, gseaParam) {
   
   # Fisher:
   if(method == "fisher") {
      geneSetStatistic <- 2*(sum(-1*log(selectedStats))) 
      
   # Stouffer & Reporter
   } else if(method %in% c("stouffer","reporter")) {
      geneSetStatistic <- sum(qnorm(selectedStats,lower.tail=FALSE)) / sqrt(length(selectedStats))
      
   # Tail strength:
   } else if(method == "tailStrength") {
      m <- length(selectedStats)
      geneSetStatistic <- (1/m)*sum(1-sort(selectedStats)*(m+1)/(1:m))
      
   # Wilcoxon rank sum:
   } else if(method == "wilcoxon_less") {
      geneSetStatistic <- unlist(wilcox.test(selectedStats, statistics, alternative="less", conf.int=FALSE)[c(1,3)])
   } else if(method == "wilcoxon_greater") {
      geneSetStatistic <- unlist(wilcox.test(selectedStats, statistics, alternative="greater", conf.int=FALSE)[c(1,3)])
   } else if(method == "wilcoxon_two.sided") {
      geneSetStatistic <- unlist(wilcox.test(selectedStats, statistics, alternative="two.sided", conf.int=FALSE)[c(1,3)])
   } else if(method == "wilcoxon_fast") {
      m <- length(selectedStats)
      #geneSetStatistic <- sum(rank(c(selectedStats,statistics[statistics <= max(selectedStats)]))[1:m])-m*(m+1)/2         
      geneSetStatistic <- sum(rank(c(selectedStats,statistics))[1:m])-m*(m+1)/2
      
   # Mean:
   } else if(method == "mean") {
      geneSetStatistic <- mean(selectedStats)
      
   # Median:   
   } else if(method == "median") {
      geneSetStatistic <- median(selectedStats)
      
   # Sum:   
   } else if(method == "sum") {
      geneSetStatistic <- sum(selectedStats)
   
   # Maxmean:
   } else if(method == "maxmean") {
      m <- length(selectedStats)
      sPlus <- sum(selectedStats[selectedStats > 0])/m
      sMinus <- -sum(selectedStats[selectedStats < 0])/m
      geneSetStatistic <- max(c(sPlus,sMinus))
      
   # GSEA:
   } else if(method == "gsea") {
      S <- selectedStats # index of genes in gene set
      r <- statistics # sorted gene-level statistics
      p <- gseaParam
      
      S <- sort(S)
      
      m <- length(S)
      N <- length(r)
      NR <- (sum(abs(r[S])^p))
      rAdj <- abs(r[S])^p
      rCumSum <- cumsum(rAdj) / NR
      
      tops <- rCumSum - (S - seq_along(S)) / (N - m) # tops/local maxima in the running sum plot
      bottoms <- tops - rAdj / NR # bottoms/local minima in the running sum plot
      maxP <- max(c(tops,0), na.rm=TRUE)
      minP <- min(c(bottoms,0), na.rm=TRUE)
      

      if(maxP > -minP) {
         geneSetStatistic <- maxP
      } else {
         geneSetStatistic <- minP
      }
      
   # PAGE:
   } else if(method == "page") {
      mu <- mean(statistics)
      delta <- stats::sd(statistics)
      Sm <- mean(selectedStats)
      m <- length(selectedStats)
      geneSetStatistic <- (Sm - mu)*sqrt(m)/delta
   }
   
   return(geneSetStatistic)
}
