\name{gsa_input}
     \docType{data}
     \alias{gsa_input}
     \title{Random input data for gene set analysis}
     \description{
       This data set is completely randomly generated and contains p-values for 2000 genes, fold-changes for those genes and a gene set collection giving the connection between genes and 50 gene sets. Only attended to be used as example data for \code{\link{runGSA}}.
     }
     \usage{data(gsa_input)}
     \format{A list containing 3 elements: gsa_input$pvals and gsa_input$directions are numeric vectors, gsa_input$gsc is a two-column matrix with gene names in the first column and gene set names in the second.}

     \keyword{datasets}
