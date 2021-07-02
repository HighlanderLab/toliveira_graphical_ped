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
library(rjags)
library(coda)
library(runjags)
library(bayesplot)
library(MCMCvis)
library(dplyr)
#=======================================================================
## Functions
#=======================================================================
file.sources <- list.files(path="./R/", pattern = "*.R")
file.sources <- paste0("./R/",file.sources)
sapply(file.sources, source, .GlobalEnv)
#=======================================================================
# Example 1 - pedigree + phenotypic records
#=======================================================================
example <- readRDS("./datasets/animal_sim10.RDS")
example$ID <- as.numeric(example$year)
#-----------------------------------------------------------------------
## For repeatability model
#-----------------------------------------------------------------------
## Hyperparameters - variance of additive genetic values and variance of
## residuals
sigma2a <- 4.5^2
h2 <- 0.4
sigma2e <- sigma2a*(1/h2 - 1)
#=======================================================================
##  Animal model - accounting inbreedin and variance known
#=======================================================================
AnimalModelBUGS

filename <- "model_example2.txt"
writeJags(model = AnimalModelBUGS, filename = filename)

## Prepair data and model
#-----------------------------------------------------------------------
## runJags
#-----------------------------------------------------------------------
jags_params <- c("a", "b")
jags_paramsVar <- c("a", "b", "sigma2a", "sigma2e")
nChains <- 3
nBurninSteps <- 5000
nThinSteps <- 10
nUseSteps <- 3000
nt <- 3 # number of replicates

tmp <- estimatesNoInbreeding <- data_cost <- ESS <- list()
for (ii in 1:6) {
  #---------------------------------------------------------------------
  # Getting the data
  #---------------------------------------------------------------------
  act_gen <- ii + 4 # actual generation
  jags.data <- example %>%
    filter(ID <= act_gen) %>%
    droplevels()
  #---------------------------------------------------------------------
  tmp[[ii]] <- AnimalModelBUGSData(
    ped=jags.data[, c(1:3, 6, 4:5)],  inbreeding = FALSE
  )

  ## Set precisions
  tmp[[ii]]$tau2a <- 1 / sigma2a
  tmp[[ii]]$tau2e <- 1 / sigma2e

  #=====================================================================
  # Jags
  #=====================================================================
  #---------------------------------------------------------------------
  # Initialization
  #---------------------------------------------------------------------
  initfunction <- function(chain) {
    nI <- tmp[[ii]]$nI
    nB <- tmp[[ii]]$nB
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
  #-----------------------------------------------------------------------
  time <- NA
  for (i in 1:nt){
    start_time <- Sys.time()
    runJagsOut <- run.jags(method = "parallel",
                           model = filename,
                           monitor = jags_params,
                           data = tmp[[ii]],
                           n.chains = nChains,
                           burnin = nBurninSteps,
                           sample = ceiling(nUseSteps/nChains),
                           thin = nThinSteps,
                           summarise = FALSE,
                           plots = FALSE,
                           inits = initfunction)
    end_time <- Sys.time()
    time[i] <- difftime(end_time, start_time, units='secs')
    cat(paste("\n", "Generation", ii, "fit number", i, "\n"))
  }
  (time_PedNoInb_mean <- mean(time))
  (time_PedNoInb_se <- sd(time)/sqrt(nt))

  # coda samples - MCMC
  coda_samples <- as.mcmc.list(runJagsOut)

  # parameter estimates
  estimatesNoInbreeding[[ii]] <- MCMCsummary(coda_samples, round = 4)

  #---------------------------------------------------------------------
  # Computational cost
  #---------------------------------------------------------------------
  neef <- estimatesNoInbreeding[[ii]]$n.eff
  cost <- as.numeric(time_PedNoInb_mean)/neef
  #---------------------------------------------------------------------
  # Summary
  #---------------------------------------------------------------------
  data_cost[[ii]] <- data.frame(time_PedNoInb_mean, time_PedNoInb_se,
                                neef, cost)
  data_cost[[ii]] <- data_cost[[ii]] %>%
    mutate_all(function(x) ifelse(is.infinite(x), NA, x))
  data_cost[[ii]]$Generation <- ii + 4
  data_cost[[ii]]$par <- gsub("[^a-zA-Z]", "",
                              rownames(estimatesNoInbreeding[[ii]]))

  ESS[[ii]] <- data_cost[[ii]] %>%
    group_by(Generation,  par) %>%
    summarise_all(
      list(mean = mean, median = median,  sd = sd,  max = max),
      na.rm = TRUE
    )
  #---------------------------------------------------------------------
}

data_cost_NoInb <- do.call(rbind.data.frame, data_cost)
data_cost_NoInb$model <- "Known Variance"

ESS_NoInb <- do.call(rbind.data.frame, ESS)
ESS_NoInb$model <- "Known Variance"

save(estimatesNoInbreeding, data_cost_NoInb,
     file = "./Results/Large_Example/DataCostNoInb.RData")
saveRDS(ESS_NoInb, file = "./Results/Large_Example/ESSNoInb.RDS")

#=======================================================================
##  Animal model - accounting inbreedin and variance unknown
#=======================================================================
AnimalModelBUGSVar

writeJags(model = AnimalModelBUGSVar, filename = filename)

#undebug(AnimalModelBUGSData)
tmp <- estimatesNoInbreedingVar <- data_cost <- ESS <- list()
for (ii in 1:6) {
  #---------------------------------------------------------------------
  # Getting the data
  #---------------------------------------------------------------------
  act_gen <- ii + 4 # actual generation
  jags.data <- example %>%
    filter(ID <= act_gen) %>%
    droplevels()
  #---------------------------------------------------------------------
  tmp[[ii]] <- AnimalModelBUGSData(
    ped=jags.data[, c(1:3, 6, 4:5)],  inbreeding = FALSE
  )
  #=====================================================================
  # Jags
  #=====================================================================
  #---------------------------------------------------------------------
  # Initialization
  #---------------------------------------------------------------------
  initfunction <- function(chain) {
    nI <- tmp[[ii]]$nI
    nB <- tmp[[ii]]$nB
    nc <- 3
    stopifnot(chain %in% 1:nc)
    a <- b <- tau2e <- tau2a <- list()
    for(i in 1:nc){
      a[[i]] <- c(rnorm(nI,0,2),NA)
      b[[i]] <- rnorm(nB,100,10)
      tau2e[[i]] <- rgamma(1, shape = 1)
      tau2a[[i]] <- rgamma(1, shape = 1)
    }
    .RNG.name <- rep("base::Super-Duper", 3)[chain]
    .RNG.seed = seq(1,3,1)[chain]
    return(list(a = a[[chain]], b = b[[chain]],
                tau2a = tau2a[[chain]], tau2e = tau2e[[chain]],
                .RNG.seed = .RNG.seed, .RNG.name = .RNG.name))
  }
  #-----------------------------------------------------------------------
  time <- NA
  for (i in 1:nt){
    start_time <- Sys.time()
    runJagsOut <- run.jags(method = "parallel",
                           model = filename,
                           monitor = jags_paramsVar,
                           data = tmp[[ii]],
                           n.chains = nChains,
                           burnin = nBurninSteps,
                           sample = ceiling(nUseSteps/nChains),
                           thin = nThinSteps,
                           summarise = FALSE,
                           plots = FALSE,
                           inits = initfunction)
    end_time <- Sys.time()
    time[i] <- difftime(end_time, start_time, units='secs')
    cat(paste("\n", "Generation", ii, "fit number", i, "\n"))
  }
  (time_PedNoInb_mean <- mean(time))
  (time_PedNoInb_se <- sd(time)/sqrt(nt))

  # coda samples - MCMC
  coda_samples <- as.mcmc.list(runJagsOut)

  # parameter estimates
  estimatesNoInbreedingVar[[ii]] <- MCMCsummary(coda_samples, round = 4)

  #---------------------------------------------------------------------
  # Computational cost
  #---------------------------------------------------------------------
  neef <- estimatesNoInbreedingVar[[ii]]$n.eff
  cost <- as.numeric(time_PedNoInb_mean)/neef
  #---------------------------------------------------------------------
  # Summary
  #---------------------------------------------------------------------
  data_cost[[ii]] <- data.frame(time_PedNoInb_mean, time_PedNoInb_se,
                                neef, cost)
  data_cost[[ii]] <- data_cost[[ii]] %>%
    mutate_all(function(x) ifelse(is.infinite(x), NA, x))
  data_cost[[ii]]$Generation <- ii + 4
  data_cost[[ii]]$par <- gsub("[^a-zA-Z]", "",
                              rownames(estimatesNoInbreedingVar[[ii]]))

  ESS[[ii]] <- data_cost[[ii]] %>%
    group_by(Generation,  par) %>%
    summarise_all(
      list(mean = mean, median = median,  sd = sd,  max = max),
      na.rm = TRUE
    )
  #---------------------------------------------------------------------
}

data_cost_NoInbVar <- do.call(rbind.data.frame, data_cost)
data_cost_NoInbVar$model <- "Unknown Variance"

ESS_NoInbVar <- do.call(rbind.data.frame, ESS)
ESS_NoInbVar$model <- "Unknown Variance"

save(estimatesNoInbreedingVar, data_cost_NoInbVar,
     file = "./Results/Large_Example/DataCostNoInbVar.RData")
saveRDS(ESS_NoInbVar, file = "./Results/Large_Example/ESSNoInbVar.RDS")
