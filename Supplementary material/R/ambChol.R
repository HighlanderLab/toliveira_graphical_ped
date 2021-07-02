#=======================================================================
# Animal Models BUGS - JAGS                                            #
# Contributors: Gregor Gorjanc,  Thiago P. Oliveira and Ivan Pocrnic   #
# Written by Gregor Gorjanc                                            #
# License: GNU General Public License version 2 (June, 1991) or later  #
# Last Update: 04 Nov 2020                                             #
#=======================================================================

#=======================================================================
# Cholesky decomposition of relationship matrix "in prior"
#=======================================================================
## Animal model in BUGS
##  - Cholesky decomposition of relationship matrix "in prior"
## Needs:
##  - Pedigree: id, fid, mid, wsqr
##  - Data: y, x, idy
##  - Constants: nI, nU, nB nY, tau2a, tau2e
## Parameters: b, u (a)

AnimalModelBUGSChol <- function()
{
  ## Additive genetic values
  for(k in 1:nI) {
    u[k] ~ dnorm(0, tau2a)
    a[id[k]] <- u[k] * wsqr[id[k]] + 0.5 * (a[fid[k]] + a[mid[k]])
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
# Cholesky decomposition of relationship matrix "in prior"
#=======================================================================
## Animal model in BUGS
##  - Cholesky decomposition of relationship matrix "in prior"
## Needs:
##  - Pedigree: id, fid, mid, wsqr
##  - Data: y, x, idy
##  - Constants: nI, nU, nB, nY, tau2a, tau2e
## Parameters: b, u (a)

AnimalModelBUGSChol2 <- function()
{
  sigmaa <- pow(tau2a, -2)

  ## Additive genetic values
  for(k in 1:nI) {
    u[k] ~ dnorm(0, 1)
    a[id[k]] <- u[k] * wsqr[id[k]] * sigmaa + 0.5 * (a[fid[k]] + a[mid[k]])
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
# Cholesky decomposition of relationship matrix "in prior" +
# variance components
#=======================================================================
#-----------------------------------------------------------------------
## Animal model in BUGS + variances
#-----------------------------------------------------------------------
##  - Cholesky decomposition of relationship matrix "in prior"
## Needs:
##  - Pedigree: id, fid, mid, wsqr
##  - Data: y, x, idy
##  - Constants: nI, nU, nB, nY
## Parameters: b, u (a), tau2a (sigma2a), tau2e (sigma2e)

AnimalModelBUGSCholVar <- function()
{
  ## Variance priors
  tau2e ~ dgamma(0.001, 0.001)
  tau2a ~ dgamma(0.001, 0.001)
  sigma2e <- 1 / tau2e
  sigma2a <- 1 / tau2a

  ## Additive genetic values
  for(k in 1:nI) {
    u[k] ~ dnorm(0, tau2a)
    a[id[k]] <- u[k] * wsqr[id[k]] + 0.5 * (a[fid[k]] + a[mid[k]])
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
#-----------------------------------------------------------------------
## Animal model in BUGS + variances
#-----------------------------------------------------------------------
##  - Cholesky decomposition of relationship matrix "in prior"
## Needs:
##  - Pedigree: id, fid, mid, wsqr
##  - Data: y, x, x, idy
##  - Constants: nI, nU, nB, nY
## Parameters: b, u (a), tau2a (sigma2a), tau2e (sigma2e)

AnimalModelBUGSChol2Var <- function()
{
  ## Variance priors
  tau2e ~ dgamma(0.001, 0.001)
  tau2a ~ dgamma(0.001, 0.001)
  sigma2e <- 1 / tau2e
  sigma2a <- 1 / tau2a
  sigmaa <- pow(tau2a, -2)

  ## Additive genetic values
  for(k in 1:nI) {
    u[k] ~ dnorm(0, 1)
    a[id[k]] <- u[k] * wsqr[id[k]] * sigmaa + 0.5 * (a[fid[k]] + a[mid[k]])
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
