\name{consensusScores}
\alias{consensusScores}
\title{
   Top consensus gene sets and boxplot
}
\description{
   Calculates the consensus scores for the gene sets using multiple gene set analysis methods (with \code{runGSA()}). Optionally also produces a boxplot to visualize the results.
}
\usage{
   consensusScores(resList, class, direction, n=50, adjusted=FALSE, method="median", 
                   plot=TRUE, cexLabel=0.8, cexLegend=1, showLegend=TRUE, rowNames="names", 
                   logScale=FALSE, main)  
}
\arguments{
   \item{resList}{
      a list where each element is an object of class \code{GSAres}, as returned by the \code{runGSA} function.
   }
   \item{class}{
      a character string determining the p-values of which directionality class that should be used as significance information for the plot. Can be one of  \code{"distinct"}, \code{"mixed"}, \code{"non"}.
   }
   \item{direction}{
      a character string giving the direction of regulation, can be either \code{"up"} or \code{"down"}.
   }
   \item{n}{
      consensus rank cutoff. All gene sets with consensus rank (see details below) \code{<=n} will be included in the plot. Defaults to 50.
   }
   \item{adjusted}{
      a logical, whether to use adjusted p-values or not. Note that if \code{runGSA} was run with the argument \code{adjMethod="none"}, the adjusted p-values will be equal to the original p-values.
   }
   \item{method}{
      a character string selecting the method, either "mean", "median", "Borda" or "Copeland".
   }
   \item{plot}{
      a logical, whether or not to draw the boxplot.
   }
   \item{cexLabel}{
      the x- and y-axis label sizes.
   }
   \item{cexLegend}{
      the legend text size.
   }
   \item{showLegend}{
      a logical, whether or not to show the legend and the indivual method ranks as points in the plot.
   }
   \item{rowNames}{
      a character string determining which rownames to use, set to either \code{"ranks"} for the consensus rank, \code{"names"} for the gene set names, or \code{"none"} to omit rownames.
   }
   \item{logScale}{
      a logical, whether or not to use log-scale for the x-axis.
   }
   \item{main}{
      a character vector giving an alternative title of the plot.
   }
}
\details{
   Based on the results given by the elements of \code{resList}, preferably representing similar runs with \code{\link{runGSA}} but with different methods, this function ranks the gene sets for each \code{GSAres} object, based on the selected directionality class. Next, the median rank for each gene set is taken as a score for top-ranking gene sets. The highest scoring gene-sets (with consensus rank, i.e. \code{rank(rankScore,ties.method="min")}, smaller or equal to \code{n}) are selected and depicted in a boxplot, showing the distribution of individual ranks (shown as colored points), as well as the median rank (shown as a red line). As an alternative of using the median rank as consensus score, it is possible to choose the mean or using the Borda or Copeland method, through the \code{method} argument.
   
   All elements of \code{resList} have to be objects containing results for the same number of gene-sets. The ranking procedure handles ties by giving them their minimum rank.

}
\value{
   A list containing a matrix of the ranks for the top \code{n} gene sets, given by each run, as well as the corresponding matrix of p-values, given by each run.
   
}
\author{
   Leif Varemo \email{piano.rpkg@gmail.com} and
   Intawat Nookaew \email{piano.rpkg@gmail.com}
}
\seealso{
   \pkg{\link{piano}}, \code{\link{runGSA}}
}
\examples{
   # Load some example GSA results:
   data(gsa_results)
      
   # Consensus scores for the top 50 gene sets (in the non-directional class):
   cs <- consensusScores(resList=gsa_results,class="non")
   
   # Access the ranks given to gene set s7 by each individual method:
   cs$rankMat["s7",]

}
