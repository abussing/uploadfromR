<p align="center">
  <img src="./images/scdeco_logo.svg" alt="scDECO logo" width="400">
</p>

<p align="center">
  <align="center"><ins>s</ins>ingle-<ins>c</ins>ell <ins>d</ins>iff<ins>e</ins>rential <ins>co</ins>-expression
</p>

<div align="center">
  <a href="https://doi.org/10.1111/biom.13701">
    <img src="https://img.shields.io/badge/DOI-doi.org%2F10.1111%2Fbiom.13701-garnet?color=%23840028" alt="DOI">
  </a>
  <a href="https://doi.org/10.1111/biom.13457">
    <img src="https://img.shields.io/badge/DOI-doi.org%2F10.1111%2Fbiom.13457-garnet?color=%23840028" alt="DOI">
  </a>
</div>





## Description

scDECO is an R package for estimating differential co-expression in single-cell RNA-seq data. 

Differential co-expression refers to correlation between two random variables $X_1, X_2$ changing with respect to the value of some covariate(s). 

<img src="images/dynamic_corr_plot.svg" alt="differential coexpression" width="700">


The package contains implementations for two different Bayesian models:
1. `scDECO.cop`: fits a bivariate Gaussian copula model with flexible, covariate-dependent, optionally zero-inflated marginals
2. `scDECO.pg`: fits a zero-inflated bivariate Poisson-Gamma model with correlation imparted through a latent bivariate normal variable


## Installation

```{r}
devtools::install_github("YenYiHo-Lab/scDECO")
library(scDECO)
```

## Usage

### `scdeco.cop`

We will the illustrate `scDECO.cop` by simulating and fitting ZINB, ZIGA data.

```{r}
n <- 2500

# simulate data
y.use <- scdeco.sim.cop(marginals=c("ZINB", "ZIGA"), x=rnorm(n),
                    eta1.true=c(-2, 0.8), eta2.true=c(-2, 0.8),
                    beta1.true=c(1, 0.5), beta2.true=c(1, 1),
                    alpha1.true=7, alpha2.true=3,
                    tau.true=c(-0.2, .3), w=runif(n,-1,1))
```
Parameters:
* `marginals`: The two marginals. Options are NB, ZINB, GA, ZIGA, Beta, ZIBEta
* `x`: The vector (or matrix) containing the covariate values to be regressed for mean and rho parameters.
* `eta1.true`: The coefficients of the 1st marginal's zero-inflation parameter. Link function is logit link.
* `eta2.true`: The coefficients of the 2nd marginal's zero-inflation parameter. Link function is logit link.
* `beta1.true`: The coefficients of the 1st marginal's mean parameter. Link function is log link.
* `beta2.true`: The coefficients of the 2nd marginal's mean parameter. Link function is log link.
* `alpha1.true`: The coefficient of the 1st marginal's second parameter. Link function is log link.
* `alpha2.true`: The coefficient of the 2nd marginal's second parameter. Link function is log link.
* `tau.true`: The coefficient of the correlation parameter. Link function is atanh(x/2).
* `w`: A vector (or matrix) containing the covariate values to be regressed for zero-inflation parameters.

```{r}
# fit the model
mcmc.out <- scdeco.cop(y=y.use, x=x.use, marginals=marginals.use, w=w.use,
                     n.mcmc=5000, burn=1000, thin=10)
```
Parameters:
* y`: 2-column matrix with the dependent variable observations.
* `n.mcmc`: The number of MCMC iterations to run.
* `burn`: The number of MCMC iterations to burn from the beginning of the chain.
* `thin`: The number of MCMC iterations to thin.

```{r
# extract estimates and confidence intervals
lowerupper <- t(apply(mcmc.out, 2, quantile, c(0.025, 0.5, 0.975)))
estmat <- cbind(lowerupper[,1],
                c(eta1.use, eta2.use, beta1.use, beta2.use, alpha1.use, alpha2.use, tau.use),
                lowerupper[,c(2,3)])
colnames(estmat) <- c("lower", "trueval", "estval", "upper")
estmat
```

And we will illustrate `scDECO.pg` by simulating correlated poisson-gamma data.

```{r}
n <- 2500

phi.use <- c(4, 4, 1/7) # over-dispersion parameters for the two marginals and the covariate
mu.use <- c(15, 15, 7) # mean parameters for the two marginals and the covariate
b.use <- c(-3, 0.1) # zero-inflation parameters
tau.use <- c(-2, 0.4) # correlation parameters

# simulate the data
simdat <- scdeco.sim.pg(N=n, b0=b.use[1], b1=b.use[2],
                        phi1=phi.use[1], phi2=phi.use[2], phi3=phi.use[3],
                        mu1=mu.use[1], mu2=mu.use[2], mu3=mu.use[3],
                        tau0=tau.use[1], tau1=tau.use[2])

# fit the model
mcmc.out <- scdeco.pg(dat=simdat,
                      b0=b.use[1], b1=b.use[2],
                      adapt_iter=500,
                      update_iter=500,
                      coda_iter=5000,
                      coda_thin=10,
                      coda_burnin=1000)

# extract the estimates and confidence intervals
estmat <- cbind(mcmc.out$quantiles[,1],
                c(1/phi.use, mu.use, tau.use),
                mcmc.out$quantiles[,c(3,5)])
colnames(estmat) <- c("lower", "true", "est", "upper")
estmat

```

## 













