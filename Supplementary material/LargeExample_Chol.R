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
AnimalModelBUGSChol

filename <- "model_example2.txt"
writeJags(model = AnimalModelBUGSChol, filename = filename)

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

tmp <- estimatesChol <- data_cost <- ESS <- list()
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
    ped=jags.data[, c(1:3, 6, 4:5)],  Chol = TRUE
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
  initfunctionChol <- function(chain) {
    nI <- tmp[[ii]]$nI
    nB <- tmp[[ii]]$nB
    nc <- 3
    stopifnot(chain %in% 1:nc)
    u <- b <- list()
    for(i in 1:nc){
      u[[i]] <- rnorm(nI,0,2)
      b[[i]] <- rnorm(nB,0,1)
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
                           inits = initfunctionChol)
    end_time <- Sys.time()
    time[i] <- difftime(end_time, start_time, units='secs')
    cat(paste("\n", "Generation", ii, "fit number", i, "\n"))
  }
  (time_PedChol_mean <- mean(time))
  (time_PedChol_se <- sd(time)/sqrt(nt))

  # coda samples - MCMC
  coda_samples <- as.mcmc.list(runJagsOut)

  # parameter estimates
  estimatesChol[[ii]] <- MCMCsummary(coda_samples, round = 4)

  #---------------------------------------------------------------------
  # Computational cost
  #---------------------------------------------------------------------
  neef <- estimatesChol[[ii]]$n.eff
  cost <- as.numeric(time_PedChol_mean)/neef
  #---------------------------------------------------------------------
  # Summary
  #---------------------------------------------------------------------
  data_cost[[ii]] <- data.frame(time_PedChol_mean, time_PedChol_se,
                                neef, cost)
  data_cost[[ii]] <- data_cost[[ii]] %>%
    mutate_all(function(x) ifelse(is.infinite(x), NA, x))
  data_cost[[ii]]$Generation <- ii + 4
  data_cost[[ii]]$par <- gsub("[^a-zA-Z]", "",
                              rownames(estimatesChol[[ii]]))

  ESS[[ii]] <- data_cost[[ii]] %>%
    group_by(Generation,  par) %>%
    summarise_all(
      list(mean = mean, median = median,  sd = sd,  max = max),
      na.rm = TRUE
    )
  #---------------------------------------------------------------------
}

data_cost_Chol <- do.call(rbind.data.frame, data_cost)
data_cost_Chol$model <- "Known Variance"

ESS_Chol <- do.call(rbind.data.frame, ESS)
ESS_Chol$model <- "Known Variance"

save(estimatesChol, data_cost_Chol,
     file = "./Results/Large_Example/DataCostChol.RData")
saveRDS(ESS_Chol, file = "./Results/Large_Example/ESSChol.RDS")

#=======================================================================
##  Animal model - accounting inbreedin and variance unknown
#=======================================================================
AnimalModelBUGSVar

writeJags(model = AnimalModelBUGSCholVar, filename = filename)

#undebug(AnimalModelBUGSData)
tmp <- estimatesCholVar <- data_cost <- ESS <- list()
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
    ped=jags.data[, c(1:3, 6, 4:5)],  Chol = TRUE
  )
  #=====================================================================
  # Jags
  #=====================================================================
  #---------------------------------------------------------------------
  # Initialization
  #---------------------------------------------------------------------
  initfunctionChol <- function(chain) {
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
                           inits = initfunctionChol)
    end_time <- Sys.time()
    time[i] <- difftime(end_time, start_time, units='secs')
    cat(paste("\n", "Generation", ii, "fit number", i, "\n"))
  }
  (time_PedChol_mean <- mean(time))
  (time_PedChol_se <- sd(time)/sqrt(nt))

  # coda samples - MCMC
  coda_samples <- as.mcmc.list(runJagsOut)

  # parameter estimates
  estimatesCholVar[[ii]] <- MCMCsummary(coda_samples, round = 4)

  #---------------------------------------------------------------------
  # Computational cost
  #---------------------------------------------------------------------
  neef <- estimatesCholVar[[ii]]$n.eff
  cost <- as.numeric(time_PedChol_mean)/neef
  #---------------------------------------------------------------------
  # Summary
  #---------------------------------------------------------------------
  data_cost[[ii]] <- data.frame(time_PedChol_mean, time_PedChol_se,
                                neef, cost)
  data_cost[[ii]] <- data_cost[[ii]] %>%
    mutate_all(function(x) ifelse(is.infinite(x), NA, x))
  data_cost[[ii]]$Generation <- ii + 4
  data_cost[[ii]]$par <- gsub("[^a-zA-Z]", "",
                              rownames(estimatesCholVar[[ii]]))

  ESS[[ii]] <- data_cost[[ii]] %>%
    group_by(Generation,  par) %>%
    summarise_all(
      list(mean = mean, median = median,  sd = sd,  max = max),
      na.rm = TRUE
    )
  #---------------------------------------------------------------------
}

data_cost_CholVar <- do.call(rbind.data.frame, data_cost)
data_cost_CholVar$model <- "Unknown Variance"

ESS_CholVar <- do.call(rbind.data.frame, ESS)
ESS_CholVar$model <- "Unknown Variance"

save(estimatesCholVar, data_cost_CholVar,
     file = "./Results/Large_Example/DataCostCholVar.RData")
saveRDS(ESS_CholVar, file = "./Results/Large_Example/ESSCholVar.RDS")
#-----------------------------------------------------------------------
# tests
#-----------------------------------------------------------------------
MCMCtrace(coda_samples[, c(5002:5009)],
          params = c(paste0("b", seq(1, 6, 1)), "sigma2a", "sigma2e"),
          ISB = FALSE,
          pdf = TRUE)
gelman.diag(coda_samples[, c(5002:5009)])

pdf("FigureChol",  width = 5,  height = 3)
autocorr.plot(coda_samples[[2]][, c(5008:5009)])
dev.off()


#=======================================================================
# Nimble
#=======================================================================
#-----------------------------------------------------------------------
# Packages
library(parallel)
#-----------------------------------------------------------------------
# Variance Known
#-----------------------------------------------------------------------
# Function
#-----------------------------------------------------------------------
run_MCMC <- function(seed = sample(1:1e5, 1),  data) {
  library(nimble)
  sigma2a <- 4.5^2
  h2 <- 0.4
  sigma2e <- sigma2a*(1/h2 - 1)
  jags_params <- c("a", "b")
  nBurninSteps <- 5000
  nThinSteps <- 10
  nUseSteps <- 1000
  #---------------------------------------------------------------------
  # Code
  #---------------------------------------------------------------------
  AnimalModelChol_Nimble <- nimbleCode({
  ## Variance priors
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
      mu[i] <- inprod(x[i, 1:nB], b[1:nB]) + a[idy[i]]
    }
  })
  #---------------------------------------------------------------------
  # Model
  #---------------------------------------------------------------------
  myModel <- nimbleModel(
    code = AnimalModelChol_Nimble,
    constants = list(nI = data$nI,
                     nY = data$nY,
                     nU = data$nU,
                     nB = data$nB,
                     id = data$id,
                     idy = data$idy,
                     fid = data$fid,
                     mid = data$mid,
                     wsqr = data$wsqr,
                     tau2a = 1/sigma2a,
                     tau2e = 1/sigma2e
                     ),
    data = list(y = data$y,
                x = data$x),
    inits = list(u = rnorm(data$nI,0,2),
                 b = rnorm(data$nB, 0, 1))
  )
  #---------------------------------------------------------------------
  # Configuring model
  #---------------------------------------------------------------------
  compModel <- compileNimble(myModel)
  confModel <- configureMCMC(
    compModel,
    monitors = jags_params
  )
  #---------------------------------------------------------------------
  # Building model
  #---------------------------------------------------------------------
  myMCMC <- buildMCMC(confModel)
  #---------------------------------------------------------------------
  # Compile MCMC
  #---------------------------------------------------------------------
  compMCMC <- compileNimble(myMCMC)
  #---------------------------------------------------------------------
  # Running MCMC
  #---------------------------------------------------------------------
  results <- runMCMC(compMCMC,
                     nburnin = nBurninSteps,
                     niter = nUseSteps * nThinSteps,
                     thin = nThinSteps,
                     setSeed = seed,
                     samplesAsCodaMCMC = TRUE)
  return(results)
}

#=======================================================================
# Paralelization
#=======================================================================
nt <- 3
tmp <- estimatesCholNimble <- data_cost <- ESS <- list()
for (ii in 1:6) {
  #---------------------------------------------------------------------
  # Getting the data
  #---------------------------------------------------------------------
  act_gen <- ii + 4 # actual generation
  jags.data <- example %>%
    filter(ID <= act_gen) %>%
    droplevels()
  #---------------------------------------------------------------------
  tmp <- AnimalModelBUGSData(
    ped=jags.data[, c(1:3, 6, 4:5)],  Chol = TRUE
  )
  #---------------------------------------------------------------------
  time <- NA
  for (i in 1:nt){
    start_time <- Sys.time()
    #-------------------------------------------------------------------
    # Model
    #-------------------------------------------------------------------
    this_cluster <- makeCluster(3)
    chains <- parLapply(cl = this_cluster, X = 1:3,
                        fun = run_MCMC,
                        data = tmp)
    stopCluster(this_cluster)
    #-------------------------------------------------------------------
    end_time <- Sys.time()
    time[i] <- difftime(end_time, start_time, units='secs')
    cat(paste("\n", "Generation", ii, "fit number", i, "\n"))
  }
  (time_PedChol_mean <- mean(time))
  (time_PedChol_se <- sd(time)/sqrt(nt))

  # coda samples - MCMC
  coda_samples <- as.mcmc.list(chains)

  # parameter estimates
  estimatesCholNimble[[ii]] <- MCMCsummary(coda_samples, round = 4)

  #---------------------------------------------------------------------
  # Computational cost
  #---------------------------------------------------------------------
  neef <- estimatesCholNimble[[ii]]$n.eff
  cost <- as.numeric(time_PedChol_mean)/neef
  #---------------------------------------------------------------------
  # Summary
  #---------------------------------------------------------------------
  data_cost[[ii]] <- data.frame(time_PedChol_mean, time_PedChol_se,
                                neef, cost)
  data_cost[[ii]] <- data_cost[[ii]] %>%
    mutate_all(function(x) ifelse(is.infinite(x), NA, x))
  data_cost[[ii]]$Generation <- ii + 4
  data_cost[[ii]]$par <- gsub("[^a-zA-Z]", "",
                              rownames(estimatesCholNimble[[ii]]))

  ESS[[ii]] <- data_cost[[ii]] %>%
    group_by(Generation,  par) %>%
    summarise_all(
      list(mean = mean, median = median,  sd = sd,  max = max),
      na.rm = TRUE
    )
  #---------------------------------------------------------------------
}

data_cost_CholNimble <- do.call(rbind.data.frame, data_cost)
data_cost_CholNimble$model <- "Known Variance"

ESS_CholNimble <- do.call(rbind.data.frame, ESS)
ESS_CholNimble$model <- "Known Variance"

save(estimatesCholNimble, data_cost_CholNimble,
     file = "./Results/Large_Example/DataCostCholNimble.RData")
saveRDS(ESS_CholNimble, file = "./Results/Large_Example/ESSCholNimble.RDS")

#-----------------------------------------------------------------------
# Variance Unknown
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
# Function
#-----------------------------------------------------------------------
run_MCMC <- function(seed = sample(1:1e5, 1),  data) {
  library(nimble)
  nBurninSteps <- 5000
  nThinSteps <- 10
  nUseSteps <- 1000
  jags_paramsVar <- c("a", "b", "sigma2a", "sigma2e")
  #---------------------------------------------------------------------
  # Code
  #---------------------------------------------------------------------
  AnimalModelChol_Nimble <- nimbleCode({
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
      mu[i] <- inprod(x[i, 1:nB], b[1:nB]) + a[idy[i]]
    }
  })
  #---------------------------------------------------------------------
  # Model
  #---------------------------------------------------------------------
  myModel <- nimbleModel(
    code = AnimalModelChol_Nimble,
    constants = list(nI = data$nI,
                     nY = data$nY,
                     nU = data$nU,
                     nB = data$nB,
                     id = data$id,
                     idy = data$idy,
                     fid = data$fid,
                     mid = data$mid,
                     wsqr = data$wsqr),
    data = list(y = data$y,
                x = data$x),
    inits = list(u = rnorm(data$nI,0,2),
                 b = rnorm(data$nB, 0, 1),
                 tau2e = rgamma(1, shape = 1),
                 tau2a = rgamma(1, shape = 1)
                 )
  )
  #---------------------------------------------------------------------
  # Configuring model
  #---------------------------------------------------------------------
  compModel <- compileNimble(myModel)
  confModel <- configureMCMC(
    compModel,
    monitors = jags_paramsVar
  )
  #---------------------------------------------------------------------
  # Building model
  #---------------------------------------------------------------------
  myMCMC <- buildMCMC(confModel)
  #---------------------------------------------------------------------
  # Compile MCMC
  #---------------------------------------------------------------------
  compMCMC <- compileNimble(myMCMC)
  #---------------------------------------------------------------------
  # Running MCMC
  #---------------------------------------------------------------------
  results <- runMCMC(compMCMC,
                     nburnin = nBurninSteps,
                     niter = nUseSteps * nThinSteps,
                     thin = nThinSteps,
                     setSeed = seed,
                     samplesAsCodaMCMC = TRUE)
  return(results)
}

#=======================================================================
# Paralelization
#=======================================================================
nt <- 3
tmp <- estimatesCholVarNimble <- data_cost <- ESS <- list()
for (ii in 1:6) {
  #---------------------------------------------------------------------
  # Getting the data
  #---------------------------------------------------------------------
  act_gen <- ii + 4 # actual generation
  jags.data <- example %>%
    filter(ID <= act_gen) %>%
    droplevels()
  #---------------------------------------------------------------------
  tmp <- AnimalModelBUGSData(
    ped=jags.data[, c(1:3, 6, 4:5)],  Chol = TRUE
  )
  #---------------------------------------------------------------------
  time <- NA
  for (i in 1:nt){
    start_time <- Sys.time()
    #-------------------------------------------------------------------
    # Model
    #-------------------------------------------------------------------
    this_cluster <- makeCluster(3)
    chains <- parLapply(cl = this_cluster, X = 1:3,
                        fun = run_MCMC,
                        data = tmp)
    stopCluster(this_cluster)
    #-------------------------------------------------------------------
    end_time <- Sys.time()
    time[i] <- difftime(end_time, start_time, units='secs')
    cat(paste("\n", "Generation", ii, "fit number", i, "\n"))
  }
  (time_PedChol_mean <- mean(time))
  (time_PedChol_se <- sd(time)/sqrt(nt))

  # coda samples - MCMC
  coda_samples <- as.mcmc.list(chains)

  # parameter estimates
  estimatesCholVarNimble[[ii]] <- MCMCsummary(coda_samples, round = 4)

  #---------------------------------------------------------------------
  # Computational cost
  #---------------------------------------------------------------------
  neef <- estimatesCholVarNimble[[ii]]$n.eff
  cost <- as.numeric(time_PedChol_mean)/neef
  #---------------------------------------------------------------------
  # Summary
  #---------------------------------------------------------------------
  data_cost[[ii]] <- data.frame(time_PedChol_mean, time_PedChol_se,
                                neef, cost)
  data_cost[[ii]] <- data_cost[[ii]] %>%
    mutate_all(function(x) ifelse(is.infinite(x), NA, x))
  data_cost[[ii]]$Generation <- ii + 4
  data_cost[[ii]]$par <- gsub("[^a-zA-Z]", "",
                              rownames(estimatesCholVarNimble[[ii]]))

  ESS[[ii]] <- data_cost[[ii]] %>%
    group_by(Generation,  par) %>%
    summarise_all(
      list(mean = mean, median = median,  sd = sd,  max = max),
      na.rm = TRUE
    )
  #---------------------------------------------------------------------
}

data_cost_CholVarNimble <- do.call(rbind.data.frame, data_cost)
data_cost_CholVarNimble$model <- "Unknown Variance"

ESS_CholVarNimble <- do.call(rbind.data.frame, ESS)
ESS_CholVarNimble$model <- "Unknown Variance"

save(estimatesCholVarNimble, data_cost_CholVarNimble,
     file = "./Results/Large_Example/DataCostCholVarNimble.RData")
saveRDS(ESS_CholVarNimble, file = "./Results/Large_Example/ESSCholVarNimble.RDS")
