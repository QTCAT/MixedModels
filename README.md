
# relMM

## Description:
An R package for (Generalized) Linear Mixed Models (G)LMM incorporating user 
defined 'G-site' covariance relationship matrices. Derived from 
[pedigreemm](https://cran.r-project.org/web/packages/pedigreemm/index.html).

This package was made, to get things done in my own research. Therefore it is 
not tested in a wide range of possible applications, but instead only for those 
areas related to my own work. This means, that you should use it only with care 
and if you know what you are doing!

## Installation:
```R
# install.packages("devtools")
devtools::install_github("jrklasen/relMM")
```

### Single observation per random-effect term levels
If a random-effect term has single observations per level, but can be 
distinguished from the residuals by the relationship matrix, the 
``control``-argument of ``relMM()`` can be used to allow a fit of this model.
```R
relMM(..., controle = lmerControl(check.nobs.vs.nlev = "ignore", 
                                  check.nobs.vs.nRE = "ignore"))
```

--------------------------------------------------------------------------------
[![License](http://img.shields.io/badge/license-GPL%20%28%3E=%202%29-brightgreen.svg?style=flat)](http://www.gnu.org/licenses/gpl-2.0.html)
