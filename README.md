# philr
*A phylogenetic transform for analysis for analysis of compositional microbiome data*

***
PhILR = Phylogenetic Isometric Log-Ratio Transform.
This R package provides functions for the analysis of compositional data (e.g., data representing proportions of different variables/parts). Specifically this package allows analysis of compositional data where the parts can be related through a phylogenetic tree (as is common in microbiota survey data) and makes available the Isometric Log Ratio transform built from the phylogenetic tree and utilizing a weighted reference measure. 

## Authors ##
Justin Silverman, Duke University 

## Citation ##
Comming soon

## License ##
All source code freely availale under [GPL-3 License](https://www.gnu.org/licenses/gpl-3.0.en.html). 

## Documentation ##
Each exported function in the package has been documented and we have also written an introductory vignette that is accessible by calling 
``` r
vignette('philr-intro', package='philr')
```

## Installation ##
If you don't already have the `devtools` package installed, run:
``` r
install.packages('devtools')
```

The development version is maintained on GitHub, try:
``` r 
devtools::install_github('jsilve24/philr')
```
Alternatively, download and decompress package and then run:
```r
source('http://bioconductor.org/biocLite.R')
devtools::install_local(‘x’)  # replace x with path to decompressed package
```

## Bugs/Feature requests ##
I appreciate bug reports and feature requests. Please post to the github issue tracker [here](https://github.com/jsilve24/philr/issues). 
