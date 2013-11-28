\name{gsaHeatmap}
\alias{gsaHeatmap}

\title{
   Heatmap of top significant gene sets
}
\description{
   This function selects the top scoring (most significant) gene sets for each directionality class and produces a heatmap plot of the results.
}
\usage{
   gsaHeatmap(gsaRes, cutoff=5, adjusted=FALSE, ncharLabel=25, cellnote="pvalue")
   }

\arguments{
  \item{gsaRes}{
      an object of class \code{GSAres}, as returned from \code{runGSA()}.
}
  \item{cutoff}{
      an integer n, so that the top n gene sets (plus possible ties) in each directionality class will be included in the heatmap.
}
   \item{adjusted}{
      a logical, whether to use adjusted p-values or not. Note that if \code{runGSA} was run with the argument \code{adjMethod="none"}, the adjusted p-values will be equal to the original p-values.
   }
   \item{ncharLabel}{
      the number of characters to include in the row labels.
   }
   \item{cellnote}{
      a character string selecting the information to be printed inside each cell of the heatmap. Either \code{"pvalue"}, \code{"rank"}, \code{"nGenes"} or \code{"none"}. Note that the actual heatmap will always be based on the consensus scores.
   }
}
\details{
   This function selects the top significant gene sets in each directionality class and draws a heatmap of the results. It provides a quick summary alternative to the \code{\link{GSAsummaryTable}} function or the \code{\link{networkPlot}}.
}
\value{
   A list, returned invisibly, containing the matrix of p-values (adjusted or non-adjusted depending on the settings) as represented in the heatmap as well as the matrix of corresponding ranks and the matrix of number of genes in each gene set (inlcuding the subset of up and down regulated genes for the mixed directional classes).
}

\author{
   Leif Varemo \email{piano.rpkg@gmail.com} and
   Intawat Nookaew \email{piano.rpkg@gmail.com}
}

\seealso{
   \pkg{\link{piano}}, \code{\link{runGSA}}, \code{\link{GSAsummaryTable}}, \code{\link{networkPlot}}
   }
\examples{

   # Load example input data to GSA:
   data("gsa_input")
   
   # Load gene set collection:
   gsc <- loadGSC(gsa_input$gsc)
      
   # Run gene set analysis:
   gsares <- runGSA(geneLevelStats=gsa_input$pvals , directions=gsa_input$directions, 
                    gsc=gsc, nPerm=500)
                    
   # Make heatmap:
   dev.new(width=10,height=10)
   gsaHeatmap(gsares)
}


