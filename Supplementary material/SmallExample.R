#=======================================================================
# Animal Models BUGS - JAGS                                            #
# Contributors: Gregor Gorjanc,  Thiago P. Oliveira and Ivan Pocrnic   #
# Written by Gregor Gorjanc                                            #
# License: GNU General Public License version 2 (June, 1991) or later  #
# Last Update: 04 Nov 2020                                             #
#=======================================================================
#=======================================================================
# Packages
#=======================================================================
if (!require("animalModels")) {
  install.packages("./R/animalModels_0.0.1.tar.gz", repos = NULL,
                   type="source")
}
library(rjags)
library(coda)
library(runjags)
library(bayesplot)
library(MCMCvis)
library(dplyr)
library(animalModels)
#=======================================================================
## Functions
#=======================================================================
file.sources <- list.files(path="./R/", pattern = "*.R")
file.sources <- paste0("./R/",file.sources)
sapply(file.sources, source, .GlobalEnv)
#=======================================================================
# Example 1 - pedigree + phenotypic records
#=======================================================================
example <- data.frame(
  individual = c(  1,   2,   3,   4,   5,   6,   7,   8,   9),
  father = c( NA,  NA,   2,   2,   4,   2,   5,  NA,   7),
  mother = c( NA,  NA,   1,  NA,   3,   3,   6,  NA,   8),
  gender = c ( 2,   1,   2,   1,   1,   2,   1,   2,   2),
  phenotype = c( NA, 105,  98, 101, 106,  93,  NA,  NA, 109),
  phenotype2 = c (NA,  NA,  10,  11,  13,   8,  NA,  NA,  15),
  group = c( NA,   1,   1,   2,   2,   2,  NA,  NA,   1)
)
example <- example[, c("individual", "father", "mother", "phenotype",
                       "group")]
#-----------------------------------------------------------------------
## For repeatability model
#-----------------------------------------------------------------------
## Hyperparameters - variance of additive genetic values and variance of
## residuals
sigma2a <- 1/3
sigma2e <- 1
(h2 <- sigma2a / (sigma2a + sigma2e)) ## 0.25
#=======================================================================
##  Animal model
#=======================================================================
AnimalModelBUGS

filename <- "model.txt"
writeJags(model = AnimalModelBUGS, filename = filename)

## Prepair data
(tmp <- AnimalModelBUGSData(ped=example))

## Set precisions
tmp$tau2a <- 1 / sigma2a
tmp$tau2e <- 1 / sigma2e

#=======================================================================
# Jags
#=======================================================================
#-----------------------------------------------------------------------
# Initialization
#-----------------------------------------------------------------------
initfunction <- function(chain) {
  nI <- 9
  nB <- 2
  nc <- 3
  stopifnot(chain %in% 1:nc)
  a <- b <- list()
  for(i in 1:nc){
    a[[i]] <- c(rnorm(nI,0,2),NA)
    b[[i]] <- rnorm(nB,100,10)
  }
  .RNG.name <- rep("base::Super-Duper", 3)[chain]
  .RNG.seed = seq(1,3,1)[chain]
  return(list(a = a[[chain]], b = b[[chain]],
              .RNG.seed = .RNG.seed, .RNG.name = .RNG.name))
}

## runJags
jags_params <- c("a", "b")
nChains <- 3
nBurninSteps <- 30000
nThinSteps <- 15
nUseSteps <- 60000
start_time <- Sys.time()
runJagsOut <- run.jags(method = "parallel",
                       model = filename,
                       monitor = jags_params,
                       data = tmp,
                       n.chains = nChains,
                       burnin = nBurninSteps,
                       sample = ceiling(nUseSteps/nChains),
                       thin = nThinSteps,
                       summarise = FALSE,
                       plots = FALSE,
                       inits = initfunction)
end_time <- Sys.time()
time_Inbreeding <- end_time - start_time

# coda samples - MCMC
coda_samples <- as.mcmc.list(runJagsOut)

# parameter estimates
estimatesInbreeding <- MCMCsummary(coda_samples, round = 4)
estimatesInbreeding

#-----------------------------------------------------------------------
# Posterior
#-----------------------------------------------------------------------
posterior <- as.matrix(coda_samples)

#-----------------------------------------------------------------------
# General Diagnostic
#-----------------------------------------------------------------------
# parameter a
mcmc_areas(posterior[, 1:9],  prob = 0.95)

# Parameter b
mcmc_areas(posterior[, 11:12],  prob = 0.95)

# trace plot
mcmc_trace(posterior)

# Autocorrelation
autocorr.plot(posterior)

# Gelman-Rubin convergence diagnostic
gelman.diag(coda_samples[, c(1:9, 11:12)])

#=======================================================================
# Animal Model ignoring inbreeding
#=======================================================================
AnimalModelBUGS

filename <- "model.txt"
writeJags(model = AnimalModelBUGS, filename = filename)

## Prepair data - ignoring inbreeding
(tmp <- AnimalModelBUGSData(ped=example,  inbreeding = FALSE))

## Set precisions
tmp$tau2a <- 1 / sigma2a
tmp$tau2e <- 1 / sigma2e

#-----------------------------------------------------------------------
# Jags
#-----------------------------------------------------------------------
initfunction <- function(chain) {
  nI <- 9
  nB <- 2
  nc <- 3
  stopifnot(chain %in% 1:nc)
  a <- b <- list()
  for(i in 1:nc){
    a[[i]] <- c(rnorm(nI,0,2),NA)
    b[[i]] <- rnorm(nB,100,10)
  }
  .RNG.name <- rep("base::Super-Duper", 3)[chain]
  .RNG.seed = seq(1,3,1)[chain]
  return(list(a = a[[chain]], b = b[[chain]],
              .RNG.seed = .RNG.seed, .RNG.name = .RNG.name))
}

## runJags
jags_params <- c("a", "b")
nChains <- 3
nBurninSteps <- 30000
nThinSteps <- 10
nUseSteps <- 60000
start_time <- Sys.time()
runJagsOut <- run.jags(method = "parallel",
                       model = filename,
                       monitor = jags_params,
                       data = tmp,
                       n.chains = nChains,
                       burnin = nBurninSteps,
                       sample = ceiling(nUseSteps/nChains),
                       thin = nThinSteps,
                       summarise = FALSE,
                       plots = FALSE,
                       inits = initfunction)
end_time <- Sys.time()
time_NoInbreeding <- end_time - start_time

# coda samples - MCMC
coda_samples <- as.mcmc.list(runJagsOut)

# parameter estimates
estimatesNoInbreeding <- MCMCsummary(coda_samples, round = 4)
estimatesNoInbreeding

#-----------------------------------------------------------------------
# Posterior
#-----------------------------------------------------------------------
posterior <- as.matrix(coda_samples)

#-----------------------------------------------------------------------
# General Diagnostic
#-----------------------------------------------------------------------
# parameter a
mcmc_areas(posterior[, 1:9],  prob = 0.95)

# Parameter b
mcmc_areas(posterior[, 11:12],  prob = 0.95)

# trace plot
mcmc_trace(posterior)

# Autocorrelation
autocorr.plot(posterior)

# Gelman-Rubin convergence diagnostic
gelman.diag(coda_samples[, c(1:9, 11:12)])
#-----------------------------------------------------------------------
# Camparing with and without inbreeding
#-----------------------------------------------------------------------
data_comp <- data.frame(NoInb = estimatesNoInbreeding[c(1:9), 1],
                        Inb = estimatesInbreeding[c(1:9), 1],
                        Animal = gl(n = 9,  k = 1,
                                    labels = paste0("Animal ", 1:9)))

data_comp %>%
  ggplot(aes(y = NoInb,  x = Inb,  label = Animal)) +
  geom_point(shape = 1,  alpha = 0.9) +
  geom_abline(intercept = 0,  slope = 1) +
  geom_text_repel(
    nudge_y      = 0.05,
    direction    = "x",
    vjust        = 0,
    segment.size = 0.2) +
  ylab("Estimates considering no inbreeding") +
  xlab("Estimates considering inbreeding") +
  theme_bw(base_size = 12)

data_diff <- data.frame(
  Parameter = rownames(estimatesInbreeding)[c(1:9, 11:12)],
  Difference = estimatesNoInbreeding[c(1:9,  11, 12), 1] -
    estimatesInbreeding[c(1:9,  11, 12), 1],
  Ratio = estimatesNoInbreeding[c(1:9,  11, 12), 1]/
    estimatesInbreeding[c(1:9,  11, 12), 1])

data_diff %>%
  kbl(
    caption = "Difference and ratio between parameter estimates",
    digits = 4,
  ) %>%
  kable_paper("hover", full_width = F) %>%
  column_spec(2, color = "white",
              background = spec_color(data_diff[,2], end = 0.7),
              popover = paste("am:", data_diff[,2])) %>%
  column_spec(3, color = "white",
              background = spec_color(data_diff[,2], end = 0.7),
              popover = paste("am:", data_diff[,3]))
#=======================================================================
## Example 3: Cholesky decomposition of relationship matrix "in prior"
#=======================================================================
writeJags(model = AnimalModelBUGSChol, filename = "model.txt")

# Database
(tmp <- AnimalModelBUGSData(ped=example, Chol=TRUE))

## Set precisions
tmp$tau2a <- 1 / sigma2a
tmp$tau2e <- 1 / sigma2e

#-----------------------------------------------------------------------
# model
#-----------------------------------------------------------------------
AnimalModelBUGSChol

initfunctionChol <- function(chain) {
  nI <- 9
  nB <- 2
  nc <- 3
  stopifnot(chain %in% 1:nc)
  u <- b <- list()
  for(i in 1:nc){
    u[[i]] <- rnorm(nI,0,2)
    b[[i]] <- rnorm(nB,100,10)
  }
  .RNG.name <- rep("base::Super-Duper", 3)[chain]
  .RNG.seed = seq(1,3,1)[chain]
  return(list(u = u[[chain]], b = b[[chain]],
              .RNG.seed = .RNG.seed, .RNG.name = .RNG.name))
}

jags_params <- c("a", "b")
start_time <- Sys.time()
runJagsOut <- run.jags(method = "parallel",
                       model = filename,
                       monitor = jags_params,
                       data = tmp,
                       n.chains = nChains,
                       burnin = nBurninSteps,
                       sample = ceiling(nUseSteps/nChains),
                       thin = nThinSteps,
                       summarise = FALSE,
                       plots = FALSE,
                       inits = initfunctionChol)
end_time <- Sys.time()
time_Chol <- end_time - start_time

# coda samples - MCMC
coda_samples <- as.mcmc.list(runJagsOut)

# parameter estimates
estimatesCholInbreeding <- MCMCsummary(coda_samples, round = 4)
estimatesCholInbreeding

#-----------------------------------------------------------------------
# Camparing models considering inbreeding -
# model 1: estimates using general inverse matrix
# model 2: estimates using Cholesky decomposition
#-----------------------------------------------------------------------
data_comp <- data.frame(CholInb = estimatesCholInbreeding[c(1:9), 1],
                        Inb = estimatesInbreeding[c(1:9), 1],
                        Animal = gl(n = 9,  k = 1,
                                    labels = paste0("Animal ", 1:9)))

data_comp %>%
  ggplot(aes(y = CholInb,  x = Inb,  label = Animal)) +
  geom_point(shape = 1,  alpha = 0.9) +
  geom_abline(intercept = 0,  slope = 1) +
  geom_text_repel(
    nudge_y      = 0.05,
    direction    = "x",
    vjust        = 0,
    segment.size = 0.2) +
  ylab("Estimates considering the Cholesky decomposition") +
  xlab("Estimates considering standard inverse matrix") +
  theme_bw(base_size = 12)

data_diff <- data.frame(
  Parameter = rownames(estimatesInbreeding)[c(1:9, 11:12)],
  Difference = estimatesCholInbreeding[c(1:9,  11, 12), 1] -
    estimatesInbreeding[c(1:9,  11, 12), 1],
  Ratio = estimatesCholInbreeding[c(1:9,  11, 12), 1]/
    estimatesInbreeding[c(1:9,  11, 12), 1])

data_diff %>%
  kbl(
    caption = "Difference and ratio between parameter estimates",
    digits = 4,
  ) %>%
  kable_paper("hover", full_width = F) %>%
  column_spec(2, color = "white",
              background = spec_color(data_diff[,2], end = 0.7),
              popover = paste("am:", data_diff[,2])) %>%
  column_spec(3, color = "white",
              background = spec_color(data_diff[,2], end = 0.7),
              popover = paste("am:", data_diff[,3]))
#-----------------------------------------------------------------------
# Reording pedegree information to evaluate inbreeding correctly
#-----------------------------------------------------------------------
writeJags(model = AnimalModelBUGSChol, filename = "model.txt")


## Prepair data
ord <- c(4, 6, 1, 9, 2, 5, 7, 8, 3)
(tmp <- AnimalModelBUGSData(ped=example[ord, ], Chol=TRUE,
                            reorder=TRUE))

## Set precisions
tmp$tau2a <- 1 / sigma2a
tmp$tau2e <- 1 / sigma2e
#-----------------------------------------------------------------------
# Model
#-----------------------------------------------------------------------
AnimalModelBUGSChol
jags_params <- c("a", "b")
start_time <- Sys.time()
runJagsOut <- run.jags(method = "parallel",
                       model = filename,
                       monitor = jags_params,
                       data = tmp,
                       n.chains = nChains,
                       burnin = nBurninSteps,
                       sample = ceiling(nUseSteps/nChains),
                       thin = nThinSteps,
                       summarise = FALSE,
                       plots = FALSE,
                       inits = initfunctionChol)
end_time <- Sys.time()
time_PedOrd <- end_time - start_time

# coda samples - MCMC
coda_samples <- as.mcmc.list(runJagsOut)

# parameter estimates
estimatesCholInbreedingOrd <- MCMCsummary(coda_samples, round = 4)
estimatesCholInbreedingOrd
#-----------------------------------------------------------------------
# Camparing models considering inbreeding
# model 2: estimates using Cholesky decomposition
# model 3: estimates using Cholesky decomposition ordering
#-----------------------------------------------------------------------
data_comp <-
  data.frame(CholInb = estimatesCholInbreeding[c(1:9), 1],
             CholInbOrd = estimatesCholInbreedingOrd[c(1:9), 1],
             Inb = estimatesInbreeding[c(1:9), 1],
             Animal = gl(n = 9,  k = 1,
                         labels = paste0("Animal ", 1:9)))

data_comp %>%
  ggplot(aes(y = CholInbOrd,  x = CholInb,  label = Animal)) +
  geom_point(shape = 1,  alpha = 0.9) +
  geom_abline(intercept = 0,  slope = 1) +
  geom_text_repel(
    nudge_y      = 0.05,
    direction    = "x",
    vjust        = 0,
    segment.size = 0.2) +
  ylab("Estimates considering the Cholesky decomposition with reorder") +
  xlab("Estimates considering the Cholesky decomposition") +
  theme_bw(base_size = 12)

data_diff <- data.frame(
  Parameter = rownames(estimatesCholInbreeding)[c(1:9, 11:12)],
  Difference = estimatesCholInbreedingOrd[c(1:9,  11, 12), 1] -
    estimatesCholInbreeding[c(1:9,  11, 12), 1],
  Ratio = estimatesCholInbreedingOrd[c(1:9,  11, 12), 1]/
    estimatesCholInbreeding[c(1:9,  11, 12), 1])

data_diff %>%
  kbl(
    caption = "Difference and ratio between parameter estimates",
    digits = 4,
  ) %>%
  kable_paper("hover", full_width = F) %>%
  column_spec(2, color = "white",
              background = spec_color(data_diff[,2], end = 0.7),
              popover = paste("am:", data_diff[,2])) %>%
  column_spec(3, color = "white",
              background = spec_color(data_diff[,2], end = 0.7),
              popover = paste("am:", data_diff[,3]))

#=======================================================================
### --- Animal model - ZStar via Cholesky ---
#=======================================================================
writeJags(model = AnimalModelBUGSZStar, filename = "model.txt")

example$group <- as.factor(example$group)
## Prepair data
(tmp <- AnimalModelBUGSData(ped=example, ZChol=TRUE, intercept = FALSE))

## Set precisions
tmp$tau2a <- 1 / sigma2a
tmp$tau2e <- 1 / sigma2e

## Store the right Cholesky factor and id for later use and remove it -
## this is mandatory!
U <- tmp$fact
tmp$fact <- NULL
id <- tmp$id
tmp$id <- NULL
#-----------------------------------------------------------------------
# Model
#-----------------------------------------------------------------------
AnimalModelBUGSZStar2
initfunctionZstar <- function(chain) {
  nI <- 6
  nB <- 2
  nc <- 3
  stopifnot(chain %in% 1:nc)
  u <- b <- list()
  for(i in 1:nc){
    u[[i]] <- rnorm(nI,0,2)
    b[[i]] <- rnorm(nB,100,10)
  }
  .RNG.name <- rep("base::Super-Duper", 3)[chain]
  .RNG.seed = seq(1,3,1)[chain]
  return(list(u = u[[chain]], b = b[[chain]],
              .RNG.seed = .RNG.seed, .RNG.name = .RNG.name))
}
jags_params <- c("u", "b")
start_time <- Sys.time()
runJagsOut <- run.jags(method = "parallel",
                       model = filename,
                       monitor = jags_params,
                       data = tmp,
                       n.chains = nChains,
                       burnin = nBurninSteps,
                       sample = ceiling(nUseSteps/nChains),
                       thin = nThinSteps,
                       summarise = FALSE,
                       plots = FALSE,
                       inits = initfunctionZstar)
end_time <- Sys.time()
time_Zstar <- end_time - start_time

# coda samples - MCMC
coda_samples <- as.mcmc.list(runJagsOut)

# parameter estimates
estimatesCholInbreedingZstar <- MCMCsummary(coda_samples, round = 4)
estimatesCholInbreedingZstar

#-----------------------------------------------------------------------
# Obtaining the a estimates
#-----------------------------------------------------------------------
tmp <- coda_samples[, 1:6]
tmp2 <- mcpar(tmp[[1]])
tmp2 <- mcmc(data=as.matrix(as.matrix(tmp) %*% U),
             start=tmp2[1], end=tmp2[2], thin=tmp2[3])
dimnames(tmp2)[[2]] <- paste("a[", id, "]", sep="")
Z_est <- c(NA, summary(tmp2)[[1]][1:5, 1], rep(NA, 2),
           summary(tmp2)[[1]][6, 1], estimatesCholInbreedingZstar[7:8, 1])
summary(tmp2)

#=======================================================================
# Animal model - ZStar via SVD
#=======================================================================
writeJags(model = AnimalModelBUGSZStar, filename = "model.txt")

## Prepair data
(tmp <- AnimalModelBUGSData(ped=example, ZSVD = TRUE))

## Set precisions
tmp$tau2a <- 1 / sigma2a
tmp$tau2e <- 1 / sigma2e

## Store the SVD factor and id for later use and remove it - this is
## mandatory!
svdUDV <- tmp$fact
tmp$fact <- NULL
id <- tmp$id
tmp$id <- NULL
#-----------------------------------------------------------------------
# Model
#-----------------------------------------------------------------------
AnimalModelBUGSZStar
jags_params <- c("u", "b")
start_time <- Sys.time()
runJagsOut <- run.jags(method = "parallel",
                       model = filename,
                       monitor = jags_params,
                       data = tmp,
                       n.chains = nChains,
                       burnin = nBurninSteps,
                       sample = ceiling(nUseSteps/nChains),
                       thin = nThinSteps,
                       summarise = FALSE,
                       plots = FALSE,
                       inits = initfunctionZstar)
end_time <- Sys.time()
time_ZStarSVD <- end_time - start_time

# coda samples - MCMC
coda_samples <- as.mcmc.list(runJagsOut)

# parameter estimates
estimatesZSVDInbreedingZstar <- MCMCsummary(coda_samples, round = 4)
estimatesZSVDInbreedingZstar

#-----------------------------------------------------------------------
# Obtaining the a estimates
#-----------------------------------------------------------------------
tmp <- coda_samples[, 1:6]
tmp2 <- mcpar(tmp[[1]])
tmp2 <- mcmc(data=as.matrix(t(svdUDV$u %*%
                                Diagonal(x=sqrt(svdUDV$d)) %*%
                                t(as.matrix(tmp[[1]])))),
             start=tmp2[1],
             end=tmp2[2],
             thin=tmp2[3])
dimnames(tmp2)[[2]] <- paste("a[", id, "]", sep="")
Z_est_ZSVD <- c(NA, summary(tmp2)[[1]][1:5, 1], rep(NA, 2),
           summary(tmp2)[[1]][6, 1], estimatesZSVDInbreedingZstar[7:8, 1])
summary(tmp2)

#=======================================================================
#=======================================================================
# # Final Summary
#=======================================================================
#=======================================================================
#-----------------------------------------------------------------------
# Comparing a estimates get from all models - time
#-----------------------------------------------------------------------
Z_est <- c(NA, summary(tmp2)[[1]][1:5, 1], rep(NA, 2),
           summary(tmp2)[[1]][6, 1], estimatesCholInbreedingZstar[7:8, 1])

Z_est_ZSVD <- c(NA, summary(tmp2)[[1]][1:5, 1], rep(NA, 2),
           summary(tmp2)[[1]][6, 1], estimatesInbreedingZstar[7:8, 1])

#-----------------------------------------------------------------------
m1 <- estimatesNoInbreeding[c(1:9, 11, 12), 1]
m2 <- estimatesInbreeding[c(1:9, 11, 12), 1]
m3 <- estimatesCholInbreeding[c(1:9, 11, 12), 1]
m4 <- estimatesCholInbreedingOrd[c(1:9, 11, 12), 1]
m5 <-  Z_est
m6 <- Z_est_ZSVD

data_comp <-
  data.frame(Parameters = rownames(estimatesCholInbreeding)[c(1:9, 11, 12)],
             m1, m2, m3, m4, m5, m6)

data_comp %>%
  kbl(
    caption = "Parameter estimates from different models",
    digits = 4,
  ) %>%
  kable_paper("hover", full_width = F)
#=======================================================================



#=======================================================================
# Computing cost
#=======================================================================
#-----------------------------------------------------------------------
# Time
#-----------------------------------------------------------------------
Time <- c(time_NoInbreeding, time_Inbreeding, time_Chol,
          time_PedOrd, time_Zstar, time_ZStarSVD)
Type <- c("No Inbreeding", "Inbreeding", "Inbreeding + Cholesky Decomposition",
          "Inbreeding + Cholesky Decomposition + Reording Pedigree",
          "Inbreeding + Cholesky Decomposition + Zstar",
          "Inbreeding + Zstar + SVD")
data_time <- data.frame(Model = paste0("m",1:6), Type= Type, "Comp Time" = Time)

data_time %>%
  kbl(
    caption = "Computational time for all models",
    digits = 4,
  ) %>%
  kable_paper("hover", full_width = F)

# -----------------------------------------------------------------------
# Comparing a estimates get from all models - effective number of
# simulation draws
# -----------------------------------------------------------------------
m1 <- estimatesNoInbreeding$n.eff[c(1:9,11,12)]
m2 <- estimatesInbreeding$n.eff[c(1:9,11,12)]
m3 <- estimatesCholInbreeding$n.eff[c(1:9,11,12)]
m4 <- estimatesCholInbreedingOrd$n.eff[c(1:9,11,12)]
m5 <- c(NA, estimatesCholInbreedingZstar$n.eff[1:5], rep(NA,2),
        estimatesCholInbreedingZstar$n.eff[6],
        estimatesCholInbreedingZstar$n.eff[7:8])
m6 <- c(NA, estimatesZSVDInbreedingZstar$n.eff[1:5], rep(NA,2),
        estimatesZSVDInbreedingZstar$n.eff[6],
        estimatesZSVDInbreedingZstar$n.eff[7:8])

neef <- rbind(m1, m2, m3, m4, m5,m6)/nUseSteps

Type <- c("No Inbreeding", "Inbreeding",
          "Inbreeding + Cholesky Decomposition",
          "Inbreeding + Cholesky Decomposition + Reording Pedigree",
          "Inbreeding + Cholesky Decomposition + Zstar",
          "Inbreeding + Zstar + SVD")
data_neef <- data.frame(Model = paste0("m",1:6), Type= Type, neef)
colnames(data_neef) <- c("Model","Type","a[1]","a[2]","a[3]","a[4]",
                         "a[5]","a[6]","a[7]","a[8]","a[9]", "b[1]",
                         "b[2]")

data_neef %>%
  kbl(
    caption = "Estimate of the effective sample size ratio ($\\widehat{R}_{eff}$) for each parameter in the model",
    digits = 4,
    row.names = FALSE
  ) %>%
  kable_styling(bootstrap_options = c("hover"), full_width = F) %>%
  add_header_above(c(" " = 2, "$\\widehat{R}_{eff}$" = 11))
#=======================================================================
