# BioFig

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
