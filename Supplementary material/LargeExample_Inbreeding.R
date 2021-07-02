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

tmp <- estimatesInbreeding <- data_cost <- ESS <- list()
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
    ped=jags.data[, c(1:3, 6, 4:5)]
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
  (time_PedInb_mean <- mean(time))
  (time_PedInb_se <- sd(time)/sqrt(nt))

  # coda samples - MCMC
  coda_samples <- as.mcmc.list(runJagsOut)

  # parameter estimates
  estimatesInbreeding[[ii]] <- MCMCsummary(coda_samples, round = 4)

  #---------------------------------------------------------------------
  # Computational cost
  #---------------------------------------------------------------------
  neef <- estimatesInbreeding[[ii]]$n.eff
  cost <- as.numeric(time_PedInb_mean)/neef
  #---------------------------------------------------------------------
  # Summary
  #---------------------------------------------------------------------
  data_cost[[ii]] <- data.frame(time_PedInb_mean, time_PedInb_se,
                                neef, cost)
  data_cost[[ii]] <- data_cost[[ii]] %>%
    mutate_all(function(x) ifelse(is.infinite(x), NA, x))
  data_cost[[ii]]$Generation <- ii + 4
  data_cost[[ii]]$par <- gsub("[^a-zA-Z]", "",
                              rownames(estimatesInbreeding[[ii]]))

  ESS[[ii]] <- data_cost[[ii]] %>%
    group_by(Generation,  par) %>%
    summarise_all(
      list(mean = mean, median = median,  sd = sd,  max = max),
      na.rm = TRUE
    )
  #---------------------------------------------------------------------
}

data_cost_Inb <- do.call(rbind.data.frame, data_cost)
data_cost_Inb$model <- "Known Variance"

ESS_Inb <- do.call(rbind.data.frame, ESS)
ESS_Inb$model <- "Known Variance"


save(estimatesInbreeding, data_cost_Inb,
     file = "./Results/Large_Example/DataCostInb.RData")
saveRDS(ESS_Inb, file = "./Results/Large_Example/ESSInb.RDS")

#=======================================================================
##  Animal model - accounting inbreedin and variance unknown
#=======================================================================
AnimalModelBUGSVar

writeJags(model = AnimalModelBUGSVar, filename = filename)

#undebug(AnimalModelBUGSData)
tmp <- estimatesInbreedingVar <- data_cost <- ESS <- list()
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
    ped=jags.data[, c(1:3, 6)], noX = TRUE
  )

  # Set Covariates
  tmp[[ii]]$x <- with(jags.data, model.matrix(~year + sex))
  tmp[[ii]]$nB <- ncol(tmp[[ii]]$x)

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
  (time_PedInb_mean <- mean(time))
  (time_PedInb_se <- sd(time)/sqrt(nt))

  # coda samples - MCMC
  coda_samples <- as.mcmc.list(runJagsOut)

  # parameter estimates
  estimatesInbreedingVar[[ii]] <- MCMCsummary(coda_samples, round = 4)

  #---------------------------------------------------------------------
  # Computational cost
  #---------------------------------------------------------------------
  neef <- estimatesInbreedingVar[[ii]]$n.eff
  cost <- as.numeric(time_PedInb_mean)/neef
  #---------------------------------------------------------------------
  # Summary
  #---------------------------------------------------------------------
  data_cost[[ii]] <- data.frame(time_PedInb_mean, time_PedInb_se,
                                neef, cost)
  data_cost[[ii]] <- data_cost[[ii]] %>%
    mutate_all(function(x) ifelse(is.infinite(x), NA, x))
  data_cost[[ii]]$Generation <- ii + 4
  data_cost[[ii]]$par <- gsub("[^a-zA-Z]", "",
                              rownames(estimatesInbreedingVar[[ii]]))

  ESS[[ii]] <- data_cost[[ii]] %>%
    group_by(Generation,  par) %>%
    summarise_all(
      list(mean = mean, median = median,  sd = sd,  max = max),
      na.rm = TRUE
    )
  #---------------------------------------------------------------------
}

data_cost_InbVar <- do.call(rbind.data.frame, data_cost)
data_cost_InbVar$model <- "Unknown Variance"

ESS_InbVar <- do.call(rbind.data.frame, ESS)
ESS_InbVar$model <- "Unknown Variance"

save(estimatesInbreedingVar, data_cost_InbVar,
     file = "./Results/Large_Example/DataCostInbVar.RData")
saveRDS(ESS_InbVar, file = "./Results/Large_Example/ESSInbVar.RDS")


#-----------------------------------------------------------------------
# tests
#-----------------------------------------------------------------------
MCMCtrace(coda_samples[, c(5002:5009)],
          params = c(paste0("b", seq(1, 6, 1)), "sigma2a", "sigma2e"),
          ISB = FALSE,
          pdf = TRUE)
gelman.diag(coda_samples[, c(5002:5009)])

pdf("FigureInb",  width = 5,  height = 3)
autocorr.plot(coda_samples[[1]][, c(5008:5009)])
dev.off()

#-----------------------------------------------------------------------
