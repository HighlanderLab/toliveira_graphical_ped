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
AnimalModelBUGSZStar

filename <- "model_example2.txt"
writeJags(model = AnimalModelBUGSZStar, filename = filename)

## Prepair data and model
#-----------------------------------------------------------------------
## runJags
#-----------------------------------------------------------------------
jags_params <- c("u", "b")
jags_paramsVar <- c("u", "b", "sigma2a", "sigma2e")
nChains <- 3
nBurninSteps <- 5000
nThinSteps <- 10
nUseSteps <- 3000
nt <- 3 # number of replicates

tmp <- estimatesCholZStar <- data_cost <- ESS <- U <- id <- list()
for (ii in 1:6) {
  #---------------------------------------------------------------------
  # Getting the data
  #---------------------------------------------------------------------
  act_gen <- ii + 4 # actual generation
  jags.data <- example %>%
    filter(ID <=3) %>%
    droplevels()
  #---------------------------------------------------------------------
  tmp[[ii]] <- AnimalModelBUGSData(
    ped=jags.data[, c(1:3, 6, 4:5)],  ZChol = TRUE
  )

  ## Set precisions
  tmp[[ii]]$tau2a <- 1 / sigma2a
  tmp[[ii]]$tau2e <- 1 / sigma2e

  ## Store the right Cholesky factor and id for later use and remove it -
  ## this is mandatory!
  U[[ii]] <- tmp[[ii]]$fact
  tmp[[ii]]$fact <- NULL
  id[[ii]] <- tmp[[ii]]$id
  tmp[[ii]]$id <- NULL

  #=====================================================================
  # Jags
  #=====================================================================
  #---------------------------------------------------------------------
  # Initialization
  #---------------------------------------------------------------------
  initfunctionZstar <- function(chain) {
    #removing ind with NA
    nI <- tmp[[ii]]$nI - sum(is.na(jags.data$phenotype))
    nB <- tmp[[ii]]$nB
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
                           inits = initfunctionZstar)
    end_time <- Sys.time()
    time[i] <- difftime(end_time, start_time, units='secs')
    cat(paste("\n", "Generation", ii, "fit number", i, "\n"))
  }
  (time_PedCholZStar_mean <- mean(time))
  (time_PedCholZStar_se <- sd(time)/sqrt(nt))

  # coda samples - MCMC
  coda_samples <- as.mcmc.list(runJagsOut)

  # parameter estimates
  estimatesCholZStar[[ii]] <- MCMCsummary(coda_samples, round = 4)

  #---------------------------------------------------------------------
  # Computational cost
  #---------------------------------------------------------------------
  neef <- estimatesCholZStar[[ii]]$n.eff
  cost <- as.numeric(time_PedCholZStar_mean)/neef
  #---------------------------------------------------------------------
  # Summary
  #---------------------------------------------------------------------
  data_cost[[ii]] <- data.frame(time_PedCholZStar_mean, time_PedCholZStar_se,
                                neef, cost)
  data_cost[[ii]] <- data_cost[[ii]] %>%
    mutate_all(function(x) ifelse(is.infinite(x), NA, x))
  data_cost[[ii]]$Generation <- ii + 4
  data_cost[[ii]]$par <- gsub("[^a-zA-Z]", "",
                              rownames(estimatesCholZStar[[ii]]))

  ESS[[ii]] <- data_cost[[ii]] %>%
    group_by(Generation,  par) %>%
    summarise_all(
      list(mean = mean, median = median,  sd = sd,  max = max),
      na.rm = TRUE
    )
  #---------------------------------------------------------------------
}

data_cost_CholZStar <- do.call(rbind.data.frame, data_cost)
data_cost_CholZStar$model <- "Known Variance"

ESS_CholZStar <- do.call(rbind.data.frame, ESS)
ESS_CholZStar$model <- "Known Variance"

save(estimatesCholZStar, data_cost_CholZStar, U,  id,
     file = "./Results/Large_Example/DataCostCholZStar.RData")
saveRDS(ESS_CholZStar, file = "./Results/Large_Example/ESSCholZStar.RDS")

#=======================================================================
##  Animal model - accounting inbreedin and variance unknown
#=======================================================================
AnimalModelBUGSVar

writeJags(model = AnimalModelBUGSCholZStarVar, filename = filename)

#undebug(AnimalModelBUGSData)
tmp <- estimatesCholZStarVar <- data_cost <- ESS <- list()
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
    ped=jags.data[, c(1:3, 6, 4:5)],  ZChol = TRUE
  )
  ## Store the right Cholesky factor and id for later use and remove it -
  ## this is mandatory!
  U[[ii]] <- tmp[[ii]]$fact
  tmp[[ii]]$fact <- NULL
  id[[ii]] <- tmp[[ii]]$id
  tmp[[ii]]$id <- NULL
  #=====================================================================
  # Jags
  #=====================================================================
  #---------------------------------------------------------------------
  # Initialization
  #---------------------------------------------------------------------
  initfunctionZStar <- function(chain) {
    nI <- tmp[[ii]]$nI
    nB <- tmp[[ii]]$nB
    nc <- 3
    stopifnot(chain %in% 1:nc)
    u <- b <- tau2e <- tau2a <- list()
    for(i in 1:nc){
      u[[i]] <- rnorm(nI,0,2)
      b[[i]] <- rnorm(nB,0,1)
      tau2e[[i]] <- rgamma(1, shape = 1)
      tau2a[[i]] <- rgamma(1, shape = 1)
    }
    .RNG.name <- rep("base::Super-Duper", 3)[chain]
    .RNG.seed = seq(1,3,1)[chain]
    return(list(u = u[[chain]], b = b[[chain]],
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
                           inits = initfunctionZStar)
    end_time <- Sys.time()
    time[i] <- difftime(end_time, start_time, units='secs')
    cat(paste("\n", "Generation", ii, "fit number", i, "\n"))
  }
  (time_PedCholZStar_mean <- mean(time))
  (time_PedCholZStar_se <- sd(time)/sqrt(nt))

  # coda samples - MCMC
  coda_samples <- as.mcmc.list(runJagsOut)

  # parameter estimates
  estimatesCholZStarVar[[ii]] <- MCMCsummary(coda_samples, round = 4)

  #---------------------------------------------------------------------
  # Computational cost
  #---------------------------------------------------------------------
  neef <- estimatesCholZStarVar[[ii]]$n.eff
  cost <- as.numeric(time_PedCholZStar_mean)/neef
  #---------------------------------------------------------------------
  # Summary
  #---------------------------------------------------------------------
  data_cost[[ii]] <- data.frame(time_PedCholZStar_mean, time_PedCholZStar_se,
                                neef, cost)
  data_cost[[ii]] <- data_cost[[ii]] %>%
    mutate_all(function(x) ifelse(is.infinite(x), NA, x))
  data_cost[[ii]]$Generation <- ii + 4
  data_cost[[ii]]$par <- gsub("[^a-zA-Z]", "",
                              rownames(estimatesCholZStarVar[[ii]]))

  ESS[[ii]] <- data_cost[[ii]] %>%
    group_by(Generation,  par) %>%
    summarise_all(
      list(mean = mean, median = median,  sd = sd,  max = max),
      na.rm = TRUE
    )
  #---------------------------------------------------------------------
}

data_cost_CholZStarVar <- do.call(rbind.data.frame, data_cost)
data_cost_CholZStarVar$model <- "Unknown Variance"

ESS_CholZStarVar <- do.call(rbind.data.frame, ESS)
ESS_CholZStarVar$model <- "Unknown Variance"

save(estimatesCholZStarVar, data_cost_CholZStarVar, U,  id,
     file = "./Results/Large_Example/DataCostCholZStarVar.RData")
saveRDS(ESS_CholZStarVar, file = "./Results/Large_Example/ESSCholZStarVar.RDS")
