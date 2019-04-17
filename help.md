---
title: piano
---

# Help

## Documentation
- [Manual/vignette](https://www.bioconductor.org/packages/release/bioc/vignettes/piano/inst/doc/piano-vignette.pdf)
- [Function reference](https://www.bioconductor.org/packages/release/bioc/manuals/piano/man/piano.pdf)

## Quick start

```
# Load piano and example data:
library(piano)
data("gsa_input")

# Load gene-set collection and take a look at the resulting object:
geneSets <- loadGSC(gsa_input$gsc)
geneSets

# Run gene-set analysis:
gsares <- rungGSA(gsa_input$pvals,
                  gsa_input$directions,
                  gsc = geneSets)

# Explore the results in an interactive Shiny app:
exploreGSAres(gsares)

```
