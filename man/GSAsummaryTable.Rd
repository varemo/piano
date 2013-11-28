\name{GSAsummaryTable}
\alias{GSAsummaryTable}
\title{
   Gene set analysis summary table
}
\description{
   Displays or saves a summary table of the results from \code{\link{runGSA}}.
}
\usage{
   GSAsummaryTable(gsaRes, save=FALSE, file=NULL)
}
\arguments{
   \item{gsaRes}{
      an object of class \code{GSAres}, as returned from \code{runGSA()}.
   }
   \item{save}{
      a logical, whether or not to save the table.
   }
   \item{file}{
      a character string giving the file name to save to.
   }
}
\details{
   The table is by default saved as an .xls file, if \code{file} is unused.
}
\value{
   The summary table as a data.frame (returned invisibly if \code{save=TRUE}).
   
}
\author{
   Leif Varemo \email{piano.rpkg@gmail.com} and
  Intawat Nookaew \email{piano.rpkg@gmail.com}
}
\seealso{
   \pkg{\link{piano}}, \code{\link{runGSA}}, \code{\link{networkPlot}}, \code{\link{gsaHeatmap}}
}
\examples{
   # Load example input data to GSA:
   data("gsa_input")
   
   # Load gene set collection:
   gsc <- loadGSC(gsa_input$gsc)
      
   # Run gene set analysis:
   gsares <- runGSA(geneLevelStats=gsa_input$pvals , directions=gsa_input$directions, 
                    gsc=gsc, nPerm=500)
      
   # Summary table:
   GSAsummaryTable(gsares)  
}


