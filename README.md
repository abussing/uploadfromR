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

Differential co-expression refers to correlation between two random variables $X_1, X_2$ changing with respect to the value of some covariate(s). 

<img src="images/dynamic_corr_plot.svg" alt="differential coexpression" width="700">


The package contains implementations for two different Bayesian models:
1. `scDECO.cop`: fits a Gaussian copula model with flexible, covariate-dependent, optionally zero-inflated marginals
2. `scDECO.pg`: fits a zero-inflated Poisson-Gamma model with correlation imparted through a latent bivariate normal variable


## Installation

```{r, eval=FALSE, message=FALSE, warning=FALSE}
devtools::install_github("YenYiHo-Lab/scDECO")
library(scDECO)
```

## Usage

We will the illustrate `scDECO.cop` by simulating and fitting ZINB, ZIGA data.

```{r}
n <- 2500

# simulate 
x.use = rnorm(n) # mu and rho covariate(s)
w.use = runif(n,-1,1) # ZINF covariate(s)

eta1.use = c(-2.2, 0.7) # ZINF parameters for 1st marginal
eta2.use = c(-2, 0.8) # ZINF parameters for 2nd marginal
beta1.use = c(1,0.5) # mean parameters for 1st marginal
beta2.use = c(1,1) # mean parameters for 2nd marginal
alpha1.use = 7 # alpha parameter for 1st marginal
alpha2.use = 3 # alpha parameter for 2nd marginal
tau.use = c(-0.2, .3) # rho parameters

# choose the marginals
marginals.use <- c("ZINB", "ZIGA")

# fit the model
y.use <- scdeco.sim.cop(marginals=marginals.use, x=x.use,
                    eta1.true=eta1.use, eta2.true=eta2.use,
                    beta1.true=beta1.use, beta2.true=beta2.use,
                    alpha1.true=alpha1.use, alpha2.true=alpha2.use,
                    tau.true=tau.use, w=w.use)
mcmc.out <- scdeco.cop(y=y.use, x=x.use, marginals=marginals.use, w=w.use,
                     n.mcmc=1000, burn=100, thin=5)

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

simdat <- scdeco.sim.pg(N=n, b0=b.use[1], b1=b.use[2],
                        phi1=phi.use[1], phi2=phi.use[2], phi3=phi.use[3],
                        mu1=mu.use[1], mu2=mu.use[2], mu3=mu.use[3],
                        tau0=tau.use[1], tau1=tau.use[2])

mcmc.out <- scdeco.pg(dat=simdat,
                      b0=b.use[1], b1=b.use[2],
                      adapt_iter=100,
                      update_iter=100,
                      coda_iter=1000,
                      coda_thin=5,
                      coda_burnin=100)

estmat <- cbind(mcmc.out$quantiles[,1],
                c(1/phi.use, mu.use, tau.use),
                mcmc.out$quantiles[,c(3,5)])
colnames(estmat) <- c("lower", "true", "est", "upper")
estmat

```








