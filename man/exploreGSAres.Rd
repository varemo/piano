\name{exploreGSAres}
\alias{exploreGSAres}
\title{
   Explore GSA results
}
\description{
   Explore GSA results interactively in a web browser using \code{shiny}.
}
\usage{
   exploreGSAres(gsaRes, browser=T, geneAnnot=NULL)
}
\arguments{
   \item{gsaRes}{
      an object of class \code{GSAres}, as returned from \code{runGSA()} or an object returned from \code{runGSAhyper()}.
   }
   \item{browser}{
      a logical, whether or not to open the Shiny app in a brwoser window. Set to \code{FALSE} to open an interactive window directly in RStudio. 
   }
   \item{geneAnnot}{
      a \code{data.frame}, containing gene annotation. The first column should be gene IDs matching those in \code{gsares}.
   }
}
\details{
  Additional gene-level information, e.g. alternative names or description, can be supplied via the \code{geneAnnot} argument. This information will show up in the gene table and the gene summary tabs.
}
\author{
   Leif Varemo \email{piano.rpkg@gmail.com}
}
\seealso{
   \pkg{\link{piano}}, \code{\link{runGSA}}, \code{\link{GSAheatmap}}
}
\examples{
   # Load example input data to GSA:
   data("gsa_input")
   
   # Load gene set collection:
   gsc <- loadGSC(gsa_input$gsc)
      
   # Run gene set analysis:
   gsares <- runGSA(geneLevelStats=gsa_input$pvals , directions=gsa_input$directions, 
                    gsc=gsc, nPerm=500)
      
   # Explore results:
   exploreGSAres(gsares)
}


