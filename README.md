## Single-species, size-structured, steady-state fish stock assessment model (s6model)

[![Tests](https://github.com/alko989/s6model/actions/workflows/install.yml/badge.svg)](https://github.com/alko989/s6model/actions/workflows/install.yml) 
[![codecov.io](https://codecov.io/github/alko989/s6model/coverage.svg?branch=rewrite)](https://codecov.io/github/alko989/s6model?branch=rewrite)

[![CRAN version](http://www.r-pkg.org/badges/version/s6model)](http://cran.r-project.org/package=s6model)

Data-poor stock assessment using the s6model


### Installation
`s6` depends on the Template Model Builder ([TMB](https://tmb-project.org))

Using `remotes' package:

``` 
install.packages('TMB', type = "source")
## install.packages('remotes')
library(remotes)

## Stable version
install_github("alko989/s6model/s6model")

## Development version
install_github("alko989/s6model/s6model", ref = "dev")
```

Another option is to download the whole code as a zipball and build and install using ```R CMD INSTALL```

### References
Kokkalis, A., Thygesen, U. H., Nielsen, A. and Andersen, K. H. (2015) ‘Limits to the reliability of size-based fishing status estimation for data-poor stocks’, Fisheries Research, 171, pp. 4–11. doi: [10.1016/j.fishres.2014.10.007](https://dx.doi.org/10.1016/j.fishres.2014.10.007).

Kokkalis, A., Eikeset, A. M., Thygesen, U. H., Steingrund, P. and Andersen, K. H. (2016) ‘Estimating uncertainty of data limited stock assessments’, ICES Journal of Marine Science: Journal du Conseil, In Press. doi: [10.1093/icesjms/fsw145](https://dx.doi.org/10.1093/icesjms/fsw145).
