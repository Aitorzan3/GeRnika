

<img src="inst/extdata/GeRnika.png" width="25%" />

`GeRnika` is an open-source R package that is capable of simulating,
visualizing and comparing tumor evolution data by using simple commands.
This aims at providing a tool to help researchers to easily simulate
tumor clonal data and analyze the results of their approaches for
studying the composition and the evolutionary history of tumors.

## Installation

`GeRnika` may be easily installed with the execution of the following
command:

``` r
# commands for installing the package from github
devtools::install_github("Aitorzan3/GeRnika")
```

Note that in order to install the vignettes together with the package,
it is necessary to attach the packages “rmarkdown”, “knitcitations”,
“knitr” and “ggpubr” to your namespace. If you have not installed them
yet, you can do it by using the following command:

``` r
# commands for installing the package together with its vignettes from github
install.packages(c("rmarkdown", "knitcitations", "knitr", "ggpubr"))
```

Once you have done that, you may use the following instruction to
install de vignettes of `GeRnika`:

``` r
# commands for installing the package together with its vignettes from github
devtools::install_github("Aitorzan3/GeRnika", build_vignettes = TRUE)
```

Once the package has been installed, the namespace of `GeRnika` may be
loaded by using

``` r
library(GeRnika)
```

This produces a short report about the versions of the packages that are
used by `GeRnika`, providing information about any conflict related to
its dependencies.

## Documentation

To view documentation for the version of this package installed in your
system, start R and enter:

``` r
browseVignettes("GeRnika")
```

## Design principles

Regarding the principles related to the design of `GeRnika`, this has
been implemented in order to be fundamentally **intuitive** and **easy
to use**. This is achieved by offering **accesible** methods for
simulating and analyzing tumor phylogenies by using **simple commands**,
**all in one** single package.
