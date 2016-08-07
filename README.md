
# relMM

### Description:
An R package for (Generalized) Linear Mixed Models (G)LMM incorporating user defined 'G-site' covariance relationship matrices (derived from [pedigreemm](https://cran.r-project.org/web/packages/pedigreemm/index.html)).

The package was created to get things done in my own research, estimation of variance components in quantitative genetics. It allows the use of methods written for the lme4 package, however, in combination with the here made changes they are perhaps not valid. Therefore, this combination should be used with caution.

### Installation:
```R
# install.packages("devtools")
devtools::install_github("jrklasen/relMM")
```

### Single observation per level of a random-effect term: 
If a random-effect term has single observations per level, but can be distinguished from the residuals by the relationship matrix, the ``control``-argument of ``relMM()`` can be used in order to bypass the checks.
```R
relMM(..., control = lmerControl(check.nobs.vs.nlev = "ignore", 
                                 check.nobs.vs.nRE = "ignore"))
```

--------------------------------------------------------------------------------
[![License](http://img.shields.io/badge/license-GPL%20%28%3E=%202%29-brightgreen.svg?style=flat)](http://www.gnu.org/licenses/gpl-2.0.html)
