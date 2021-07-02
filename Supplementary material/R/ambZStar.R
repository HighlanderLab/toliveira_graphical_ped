#=======================================================================
# Animal Models BUGS - JAGS                                            #
# Contributors: Gregor Gorjanc,  Thiago P. Oliveira and Ivan Pocrnic   #
# Written by Gregor Gorjanc                                            #
# License: GNU General Public License version 2 (June, 1991) or later  #
# Last Update: 04 Nov 2020                                             #
#=======================================================================

#-----------------------------------------------------------------------
## Animal model in BUGS
#-----------------------------------------------------------------------
##  - Decomposition of relationship matrix "in design matrix"
## Needs:
##  - Pedigree: ZStar
##  - Data: y, x
##  - Constants: nI, nB, nY, tau2a, tau2e
## Parameters: b, u (a)
AnimalModelBUGSZStar <- function()
{
  ## Transformed additive genetic values
  for(k in 1:nI) { u[k] ~ dnorm(0, tau2a) }

  ## Fixed Effects
  for(j in 1:nB) { b[j] ~ dnorm(0, 1.0E-6) }

  ## Phenotypes
  for(i in 1:nY) {
    y[i] ~ dnorm(mu[i], tau2e)
    mu[i] <- inprod(x[i, ], b[]) + inprod(ZStar[i, ], u[])
  }
}

#-----------------------------------------------------------------------
## Animal model in BUGS + variances
#-----------------------------------------------------------------------
##  - Decomposition of relationship matrix "in design matrix"
## Needs:
##  - Pedigree: ZStar
##  - Data: y, x
##  - Constants: nI, nB, nY
## Parameters: b, u (a), tau2a (sigma2a), tau2e (sigma2e)

AnimalModelBUGSZStarVar <- function()
{
  ## Variance priors
  tau2e ~ dgamma(0.001, 0.001)
  tau2a ~ dgamma(0.001, 0.001)
  sigma2e <- 1 / tau2e
  sigma2a <- 1 / tau2a

  ## Transformed additive genetic values
  for(k in 1:nI) { u[k] ~ dnorm(0, tau2a) }

  ## Fixed Effects
  for(j in 1:nB) { b[j] ~ dnorm(0, 1.0E-6) }

  ## Phenotypes
  for(i in 1:nY) {
    y[i] ~ dnorm(mu[i], tau2e)
    mu[i] <- inprod(x[i, ], b[]) + inprod(ZStar[i, ], u[])
  }
}
