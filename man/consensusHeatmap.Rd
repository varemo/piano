\name{consensusHeatmap}
\alias{consensusHeatmap}

\title{
   Heatmap of top consensus gene sets
}
\description{
   Based on multiple result objects from the \code{\link{runGSA}} function, this function computes the consensus scores, based on rank aggregation, for each directionality class and produces a heatmap plot of the results.
}
\usage{
   consensusHeatmap(resList, method = "median", cutoff = 5, adjusted=FALSE, ncharLabel=25)
}

\arguments{
  \item{resList}{
      a list where each element is an object of class \code{GSAres}, as returned by the \code{runGSA} function.
}
  \item{method}{
      a character string selecting the method, either "mean", "median", "Borda" or "Copeland".
}
  \item{cutoff}{
      the maximum consensus score of a gene set, in any of the directionality classes, to be included in the heatmap.
}
   \item{adjusted}{
      a logical, whether to use adjusted p-values or not. Note that if \code{runGSA} was run with the argument \code{adjMethod="none"}, the adjusted p-values will be equal to the original p-values.
   }
   \item{ncharLabel}{
      the number of characters to include in the row labels.
   }
}
\details{
   This function computes the consensus gene set scores for each directionality class based on the results (gene set p-values) listed in \code{resList}, using the \code{consensusScores} function. For each class, only the \code{GSAres} objects in \code{resList} that contain p-values for that class are used as a basis for the rank aggregation. Hence, if not all classes are covered by at least 2 \code{GSAres} objects in the list, the \code{consensusHeatmap} function will not work. The results are displayed in a heatmap showing the consensus scores.
}
\value{
   A list, returned invisibly, containing the matrix of consensus scores as represented in the heatmap as well as the matrix of corresponding median p-values.
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
   
   # Consensus heatmap:
   dev.new(width=10,height=10)
   consensusHeatmap(resList=gsa_results)
   
   # Store the output:
   dev.new(width=10,height=10)
   ch <- consensusHeatmap(resList=gsa_results)
   
   # Access the median p-values for gene set s1:
   ch$pMat["s1",]
}


