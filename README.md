<p align="center">
  <img src="./images/scdeco_logo.svg" alt="scDECO logo" width="400">
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

scDECO is an R package for modeling <ins>s</ins>ingle-<ins>c</ins>ell <ins>d</ins>iff<ins>e</ins>rential <ins>co</ins>-expression.

Differential co-expression refers to the phenomenon of correlation between two random variables $\boldsymbol{X}_1, \boldsymbol{X}_2$ being covariate-dependent. 

<img src="images/dynamic_corr_plot.svg" alt="Differential Co-expression" width="600">

The above plot demonstrates how $\text{corr}\left(\boldsymbol{X}_1, \boldsymbol{X}_2\right)$ goes from negative to positive as the value of $\boldsymbol{X}_3$ increases.





## Installation

```{r, eval=FALSE, message=FALSE, warning=FALSE}
devtools::install_github("YenYiHo-Lab/scDECO")
```









