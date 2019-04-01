---
title: piano
---

# An R/Bioconductor package for gene set analysis

## Overview of piano

Piano is used for running gene set analysis (GSA) using a selection of available methods (some of which are also implemented in separate packages or other software) and starting from different kinds of gene level statistics. The advantage with using piano is that all of these methods can be run using the same setup, resulting in minimal effort for the user when switching between methods. Additionally, piano contains a new approach which divides the gene set results into directionality classes,
giving deeper information about the underlying gene pattern. Further on, piano contains functions to perform consensus analysis, i.e. combining the results of multiple GSA runs. For details, see the publication below.

## Citation

Please cite the piano publication if you use it in your research:

Väremo, L., Nielsen, J. and Nookaew, I. (2013) Enriching the gene set analysis of genome-wide data by incorporating directionality of gene expression and combining statistical hypotheses and methods. *Nucleic Acids Research.* 41 (8), 4378-4391. [https://doi.org/10.1093/nar/gkt111](https://doi.org/10.1093/nar/gkt111)

## Contact

Email the developer at piano.rpkg (at) gmail.com.
If you encountered any problems please check if your questions are already answered on the help page, in the vignette or in the function documentation before you email for support. 

You may also file an [issue](https://github.com/varemo/piano/issues) or [pull request](https://github.com/varemo/piano/pulls) directly on GitHub.

## Acknowledgements

This work was carried out at [Systems and Synthetic Biology](https://sysbio.se) at [Chalmers University of Technology](https://chalmers.se), Gothenburg, Sweden. Funding was provided from Knut and Alice Wallenberg foundation, Chalmers foundation and Bioinformatics Infrastructure for Life Sciences (BILS). The computations were performed on resources provided by the Swedish National Infrastructure for Computing (SNIC) at C3SE.

Current development and maintenance is carried out in affiliation with the [National Bioinformatics Infrastructure Sweden (NBIS)](https://nbis.se).
