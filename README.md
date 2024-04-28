<p align="center">
  <img src="./images/scdeco_logo.svg" alt="scDECO logo" width="400">
</p>

<p align="center">
  <align="center"><ins>s</ins>ingle-<ins>c</ins>ell <ins>d</ins>iff<ins>e</ins>rential <ins>co</ins>-expression
</p>

<div align="center">
  <a href="https://doi.org/10.1111/biom.13701">
    <img src="https://img.shields.io/badge/DOI-doi.org%2F10.1111%2Fbiom.13701-blue" alt="DOI">
  </a>
  <a href="https://doi.org/10.1111/biom.13457">
    <img src="https://img.shields.io/badge/DOI-doi.org%2F10.1111%2Fbiom.13457-blue" alt="DOI">
  </a>
</div>





## Description

scDECO is an R package for estimating differential co-expression in single-cell RNA-seq data. 

The package contains implementations for two different Bayesian models:
1. `scDECO.cop`: A Gaussian copula model with flexible, covariate-dependent, optionally zero-inflated marginals. 
2. `scDECO.pg`: A zero-inflated Poisson-Gamma model with non-covariate-dependent negative binomial marginals. 




Differential co-expression refers to covariate-dependence in the correlation between two random variables $X_1, X_2$. 

For exampthe correlation between $X_1$ and $X_2$ changes from strongly negative to strongly positive as the value of covariate $X_3$ grows larger.

<img src="images/dynamic_corr_plot.svg" alt="differential coexpression" width="600">

## Installation

```{r, eval=FALSE, message=FALSE, warning=FALSE}
devtools::install_github("YenYiHo-Lab/scDECO")
```









