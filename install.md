---
title: piano
---

# Install piano

## Stable release version
This is suitable for most users. Install in R:
```
BiocManager::install("piano")
```
See also [piano on Bioconductor.org](https://www.bioconductor.org/packages/release/bioc/html/piano.html) and the [changelog](http://www.bioconductor.org/packages/release/bioc/news/piano/NEWS) for more information.

## Development version
For the latest features, install the Bioconductor development version according to [these instructions](http://bioconductor.org/developers/how-to/useDevel/). Also see [this guide](https://cran.r-project.org/web/packages/BiocManager/vignettes/BiocManager.html) on how to use `BiocManager` with different package/Bioconductor/R versions.

The cutting edge code is available in the master branch (in principle identical to the Bioconductor development version) and in feature branches on [GitHub](https://github.com/varemo/piano). Install directly from GitHub in R with the command: 
```
devtools::install_github("varemo/piano", ref="<branch-name>")
```

## Conda

Piano is available through the [Bioconda channel](http://bioconda.github.io/recipes/bioconductor-piano/README.html). Use the command below to install:

```
conda install bioconductor-piano
```
For reproducibility it is recommended to explicitly install a specific version, e.g. `conda install bioconductor-piano:1.22.0`.

## Docker

Bioconda packages have a corresponding Docker BioContainer automatically created and uploaded to [Quay.io](https://Quay.io). Pull the image for piano:
```
docker pull quay.io/biocontainers/bioconductor-piano:<tag>
```
Note: the default tag `latest` might not work. See [bioconductor-piano/tags](https://quay.io/repository/biocontainers/bioconductor-piano?tab=tags) for valid values for `<tag>`.

In addition, Bioconductor also supplies [Docker images](https://www.bioconductor.org/help/docker/) for running R and Bioconductor packages.
