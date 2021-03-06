odmgm
================

<!-- README.md is generated from README.Rmd. Please edit that file -->
[![DOI](https://zenodo.org/badge/290208715.svg)](https://zenodo.org/badge/latestdoi/290208715)

This research compendium contains `R` code used to conduct the analysis of our work \[reference to paper\]. The contents are organized as an `R` package and it is therefore possible to install the package, called `odmgm`, as follows. First run the following command

``` r
install.packages("devtools")
devtools::install_github("mlindsk/odmgm")
```

Secondly, the package suggests the `gRapHD` package, all though no longer maintained on `CRAN`. To install `gRapHD`, first make a git clone, and then call the following command

``` r
install.packages("gRapHD/gRapHD_0.2.6.tar.gz", repos = NULL, type = "source")
```

The `gRapHD` package is needed to fit a mixed interaction graph. Currently, the `R` package ecosystem lacks good maintained algorithms for the purpose of model selection in mixed graphical models. We provide the wrapper function `odmgm:::fit_mixed_graph`.

The Analysis
============

The code snippets to reproduce the results in the paper is located in the `analysis` folder.

See also
========

For outlier detection in pure categorical data, see the `R` package (on `CRAN`) [molic](https://github.com/mlindsk/molic). The `molic` package depends on the `R` package (on `CRAN`) [ess](https://github.com/mlindsk/ess) which is designed for efficient model selection in discrete decomposable models.
