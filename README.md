# GeRnika
GeRnika is a package capable of simulating tumor data, visualizing it by means of phylogenetic trees and comparing tumor phylogenies. GeRnika aims at providing a tool to help researchers to easily simulate tumor data and analyze the results of their approaches for studying the composition and the evolutionary history of tumors.

# GeRnika package

`GeRnika` may be easily installed with the execution of the following R commands:

```{r, eval = FALSE}
# commands for installing the package
devtools::install_github("Aitorzan3/GeRnika", build_vignettes = TRUE)
```
Once installed, the namespace of `GeRnika` may be loaded by attaching it:

```{r setup}
library(GeRnika)
```

This produces a short report about the versions of the packages that are used by `GeRnika`, providing information about any conflicts with previously loaded packages. `GeRnika` includes the following packages: data.tree [@2], tidyverse [@3], Diagrammer [@4], MCMCpack [@5], reshape2 [@6] and colorspace [@11]

# Design principles
Regarding the principles related to the design of `GeRnika`, this has been implemented in order to be fundamentally **intuitive** and **easy to use**. There exists different *R* packages that allow their users to solve the Cloning Deconvolution Problem (CDP) and analyze the phylogeny of tumor samples, but they do not offer many approaches for visualizing phylogenetic trees in a comprehensive way nor compare the phylogeny of different samples. 

Following the above, `GeRnika` offers **accesible** methods for simulating and analyzing tumors phylogeny by using **simple commands**, **all in one** single package.

# Acknowledgments

`GeRnika` would not have been possible without the hard work of the *Intelligent Systems Group* members, specially **Borja Calvo** and **Maitena Tellaetxe**, who implemented various methods gathered in this package and coordinated its creation. 
