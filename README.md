
# MixedModels

### Description:
An R package for (Generalized) Linear Mixed Models (G)LMM incorporating user defined 'G-site' covariance relationship matrices (derived from [pedigreemm](https://cran.r-project.org/web/packages/pedigreemm/index.html)).

The package was created in order to get things done in my own research; estimation of variance components in quantitative genetics. Any other use case is not tested and may or may not be valid.

### Installation:
```R
# install.packages("devtools")
devtools::install_github("jrklasen/MixedModels")
```

### Single observation per level of a random-effect term: 
If a random-effect term has single observations per level, but can be distinguished from the residuals by the relationship matrix, the ``control``-argument of ``mm()`` can be used in order to bypass the checks.
```R
mm(..., control = lmerControl(check.nobs.vs.nlev = "ignore", 
                              check.nobs.vs.nRE = "ignore"))
```

--------------------------------------------------------------------------------
[![License](http://img.shields.io/badge/license-GPL%20%28%3E=%202%29-brightgreen.svg?style=flat)](http://www.gnu.org/licenses/gpl-2.0.html)
