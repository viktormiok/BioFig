# tigaRcycle

The R-package `tigaRcycle` performs integrative detection of circadian signals in time series omics data, both sequencing counts RNAseq and continuous microarray data. The package offers functions for visualisation of 


Note: if you have a choice to use either Windows or Unix/Linux, opt for the latter. `tigaR` runs more efficiently under Unix/Linux than under Windows. NOTE:  when running `tigaR` you may see *** WARNINGS ***  from `INLA` (e.g. on eigenvalues, or on convergence, or even something like 18500 Aborted...). They can currently not be surpressed, because they are produced by C-code. Please ignore them. 

# Installation

The package `tigaRcycle` depends on [tigaR](https://github.com/viktormiok/tigaR), [R >= 3.0.0](https://cran.r-project.org/) and is available from GitHub. This requires the package [devtools](https://cran.r-project.org/web/packages/devtools/index.html):

``` r
devtools::install_github("viktormiok/tigaRcycle", build_vignettes = TRUE)
```

Please restart R before loading the package and its documentation:

``` r
library(tigaRcycle)
utils::help(tigaRcycle)
utils::vignette("tigaRcycle")
```

# References

Publication related to `tigaRcycle` include:

- Miok, V., Wilting, S.M., Van de Wiel, M.A., Jaspers, A., van Noort, P.I., Brakenhoff, R.H., Snijders, P.J.F., Steenbergen, R.D.M., Van Wieringen, W.N. (2014), "tigaR: integrative significance analysis of temporal differential gene expression induced by genomic abnormalities", *BMC Bioinformatics*, 15, 327. ([doi.org/10.1186/1471-2105-15-327](https://doi.org/10.1186/1471-2105-15-327)).

Please cite the relevant publication if you use `tigaRcycle`.
