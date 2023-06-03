<img src="https://github.com/viktormiok/viktormiok.me/blob/main/software/BioFig.png" align="right" height="200" width="200">

![](https://img.shields.io/badge/language-R-orange.svg) ![version](https://img.shields.io/badge/GiHub_version-1.1.0-519dd9) ![GitHub code size in bytes](https://img.shields.io/github/languages/code-size/viktormiok/tigaR) ![GitHub issues](https://img.shields.io/github/issues/viktormiok/tigaR)

![dependencies](https://img.shields.io/badge/dependencies-up%20to%20date-orange)  	![commit](https://img.shields.io/github/last-commit/viktormiok/tigaR) ![GitHub](https://img.shields.io/github/license/viktormiok/tigaR)

[![Edit with Gitpod](https://gitpod.io/button/open-in-gitpod.svg)](https://gitpod.io/#https://github.com/viktormiok/tigaR) 



# BioFig

- [Overview](#overview)
  * [Application](#application)
- [Installation](#installation)
- [Docker](#docker)
- [Data](#data)
- [Tutorials](#tutorials)
- [License](#license)
- [References](#references)

## Overview
The R-package `BioFig` performs visualisatio and quality contro of omics data. 

# Installation

The package `BioFig` depends on [BioFig](https://github.com/viktormiok/BioFig), [R >= 3.0.0](https://cran.r-project.org/) and is available from GitHub. This requires the package [devtools](https://cran.r-project.org/web/packages/devtools/index.html):

``` r
devtools::install_github("viktormiok/BioFig", build_vignettes = TRUE)
```

Please restart R before loading the package and its documentation:

``` r
library(BioFig)
utils::help(BioFig)
utils::vignette("BioFig")
```
