#=======================================================================
# Animal Models BUGS - JAGS                                            #
# Contributors: Gregor Gorjanc,  Thiago P. Oliveira and Ivan Pocrnic   #
# Written by Gregor Gorjanc                                            #
# License: GNU General Public License version 2 (June, 1991) or later  #
# Last Update: 04 Nov 2020                                             #
#=======================================================================
#=======================================================================
## Animal model in BUGS + estimation of variances components
#=======================================================================
## Needs:
##  - Pedigree: id, fid, mid, winv
##  - Data: y, x, idy
##  - Constants: nI, nU, nB, nY
## Parameters: b, a, sigma_a = 1/tau2a, sigma_e = 1/tau2e

AnimalModelBUGSVar <- function()
{
  ## Priors - precision parameters
  tau2e ~ dgamma(0.001, 0.001)
  tau2a ~ dgamma(0.001, 0.001)
  sigma2e <- pow(tau2e, -1)
  sigma2a <- pow(tau2a, -1)

  ## Additive genetic values
  for(k in 1:nI) {
    a[id[k]] ~ dnorm(pa[id[k]], Xtau2a[id[k]])
    pa[id[k]] <- 0.5 * (a[fid[k]] + a[mid[k]])
    Xtau2a[id[k]] <- winv[id[k]] * tau2a
  }
  a[nU] <- 0 # NULL (zero) holder

  ## Fixed Effects
  for(j in 1:nB) { b[j] ~ dnorm(0, 1.0E-6) }

  ## Phenotypes
  for(i in 1:nY) {
    y[i] ~ dnorm(mu[i], tau2e)
    mu[i] <- inprod(x[i, ], b[]) + a[idy[i]]
  }
}

#=======================================================================
#  Animal model in BUGS + estimation of variances components -
# prior to sigma
#=======================================================================
## Animal model in BUGS + variances
## Needs:
##  - Pedigree: id, fid, mid, winv
##  - Data: y, x, idy
##  - Constants: nI, nU, nB, nY
## Parameters: b, a, tau2a (sigmaa), tau2e (sigmae)
AnimalModelBUGSVar2 <- function()
{
 ## Variance priors
  sigmae ~ dunif(0, sigmaeu)
  sigmaa ~ dunif(0, sigmaau)
  sigma2e <- sigmae * sigmae
  sigma2a <- sigmaa * sigmaa
  tau2e <- 1 / sigma2e
  tau2a <- 1 / sigma2a

  ## Additive genetic values
  for(k in 1:nI) {
    a[id[k]] ~ dnorm(pa[id[k]], Xtau2a[id[k]])
    pa[id[k]] <- 0.5 * (a[fid[k]] + a[mid[k]])
    Xtau2a[id[k]] <- winv[id[k]] * tau2a
  }
  a[nU] <- 0 # NULL (zero) holder

  ## Fixed Effects
  for(j in 1:nB) { b[j] ~ dnorm(0, 1.0E-6) }

  ## Phenotypes
  for(i in 1:nY) {
    y[i] ~ dnorm(mu[i], tau2e)
    mu[i] <- inprod(x[i, ], b[]) + a[idy[i]]
  }
}

#=======================================================================
# Animal model in BUGS - variance terms are known a-priori
#=======================================================================
## Animal model in BUGS
## Needs:
##  - Pedigree: id, fid, mid, winv
##  - Data: y, x, idy
##  - Constants: nI, nU, nB, nY, tau2a, tau2e
## Parameters: b, a

AnimalModelBUGS <- function()
{
  ## Additive genetic values
  for(k in 1:nI) {
    a[id[k]] ~ dnorm(pa[id[k]], Xtau2a[id[k]])
    pa[id[k]] <- 0.5 * (a[fid[k]] + a[mid[k]])
    Xtau2a[id[k]] <- winv[id[k]] * tau2a
  }
  a[nU] <- 0 # NULL (zero) holder

  ## Fixed Effects
  for(j in 1:nB) { b[j] ~ dnorm(0, 1.0E-6) }

  ## Phenotypes
  for(i in 1:nY) {
    y[i] ~ dnorm(mu[i], tau2e)
    mu[i] <- inprod(x[i, ], b[]) + a[idy[i]]
  }
}
#=======================================================================
# Variant 1
#=======================================================================
AnimalModelBUGSVariant1Var <- function()
{
  ## Priors - precision parameters
  tau2e ~ dgamma(0.001, 0.001)
  tau2a ~ dgamma(0.001, 0.001)
  sigma2e <- pow(tau2e, -1)
  sigma2a <- pow(tau2a, -1)

  ## Transformed additive genetic values - variant 1
  for(k in 1:nI) {
    u[k] ~ dnorm(0, tau2a)
    a[id[k]] <- u[k] * wsqr[id[k]] +
      0.5 * (a[fid[k]] + a[mid[k]])
    }
  a[nU] <- 0 # NULL (zero) holder

  ## Location priors
  for(j in 1:nB) {
    b[j] ~ dnorm(0, 1.0E-6)
  }

  ## Phenotypes
  for(i in 1:nY) {
    y[i] ~ dnorm(mu[i], tau2e)
    mu[i] <- inprod(x[i, ], b[]) + a[idy[i]]
  }
}

AnimalModelBUGSVariant1 <- function()
{
  ## Transformed additive genetic values - variant 1
  for(k in 1:nI) {
    u[k] ~ dnorm(0, tau2a)
    a[id[k]] <- u[k] * wsqr[id[k]] +
      0.5 * (a[fid[k]] + a[mid[k]])
    }
  a[nU] <- 0 # NULL (zero) holder

  ## Location priors
  for(j in 1:nB) {
    b[j] ~ dnorm(0, 1.0E-6)
  }

  ## Phenotypes
  for(i in 1:nY) {
    y[i] ~ dnorm(mu[i], tau2e)
    mu[i] <- inprod(x[i, ], b[]) + a[idy[i]]
  }
}
#=======================================================================
# Variant 2
#=======================================================================
AnimalModelBUGSVariant2Var <- function()
{
  ## Priors - precision parameters
  tau2e ~ dgamma(0.001, 0.001)
  tau2a ~ dgamma(0.001, 0.001)
  sigma2e <- pow(tau2e, -1)
  sigma2a <- pow(tau2a, -1)

  ## Transformed additive genetic values - variant 2
  sigmaa <- pow(tau2a, -2)
  for(k in 1:nI) {
    u[k] ~ dnorm(0, 1)
    a[id[k]] <- u[k] * wsqr[id[k]] * sigmaa +
      0.5 * (a[fid[k]] + a[mid[k]])
    }
  a[nU] <- 0 # NULL (zero) holder

  ## Location priors
  for(j in 1:nB) {
    b[j] ~ dnorm(0, 1.0E-6)
  }

  ## Phenotypes
  for(i in 1:nY) {
    y[i] ~ dnorm(mu[i], tau2e)
    mu[i] <- inprod(x[i, ], b[]) + a[idy[i]]
  }
}

AnimalModelBUGSVariant2 <- function()
{
  ## Transformed additive genetic values - variant 2
  sigmaa <- pow(tau2a, -2)
  for(k in 1:nI) {
    u[k] ~ dnorm(0, 1)
    a[id[k]] <- u[k] * wsqr[id[k]] * sigmaa +
      0.5 * (a[fid[k]] + a[mid[k]])
    }
  a[nU] <- 0 # NULL (zero) holder

  ## Location priors
  for(j in 1:nB) {
    b[j] ~ dnorm(0, 1.0E-6)
  }

  ## Phenotypes
  for(i in 1:nY) {
    y[i] ~ dnorm(mu[i], tau2e)
    mu[i] <- inprod(x[i, ], b[]) + a[idy[i]]
  }
}
#=======================================================================
