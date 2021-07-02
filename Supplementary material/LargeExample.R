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
library(ggrepel)
library(dplyr)
library(kableExtra)
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

filename <- "model.txt"
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

#undebug(AnimalModelBUGSData)
tmp <- runJagsModelInb <- estimatesInbreeding <-
   data_cost <- ESS <- list()
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
    time[i] <- end_time - start_time
    cat(paste("\n", "Generation", ii, "fit number", i, "\n"))
  }
  (time_PedInb_mean <- mean(time))
  (time_PedInb_se <- sd(time)/sqrt(nt))

  # coda samples - MCMC
  coda_samples <- as.mcmc.list(runJagsOut)

  # parameter estimates
  runJagsModelInb[[ii]] <- runJagsOut
  estimatesInbreeding[[ii]] <- MCMCsummary(coda_samples, round = 4)

  #---------------------------------------------------------------------
  # Computational cost
  #---------------------------------------------------------------------
  par <- rownames(estimatesInbreeding[[ii]])
  neef <- estimatesInbreeding[[ii]]$n.eff
  cost <- as.numeric(time_PedInb_mean)/neef
  data_cost[[ii]] <- data.frame(par, time_PedInb_mean, time_PedInb_se,
                                neef, rho, cost)
  data_cost[[ii]] <- data_cost[[ii]] %>%
    mutate_all(function(x) ifelse(is.infinite(x), NA, x))
  data_cost[[ii]]$Generation <- ii + 4
  ESS[[ii]] <- data_cost[[ii]] %>%
    summarise_all(
      list(mean = mean, median = median,  sd = sd), na.rm = TRUE
    )
  #---------------------------------------------------------------------
}

data_cost_Inb <- do.call(rbind.data.frame, data_cost)
ESS_Inb <- do.call(rbind.data.frame, ESS)

save(runJagsModelInb, estimatesInbreeding, data_cost_Inb, ESS_Inb,
     file = paste0("./JAGS/runJagsOutInb",ii, ".RData"))

#=======================================================================
##  Animal model - accounting inbreedin and variance unknown
#=======================================================================
AnimalModelBUGSVar

filename <- "model.txt"
writeJags(model = AnimalModelBUGSVar, filename = filename)

#undebug(AnimalModelBUGSData)
tmp <- runJagsModelInbVar <- estimatesInbreedingVar <-
   data_costInbVar <- ESSInbVar <- list()
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
  runJagsModelInbVar[[ii]] <- runJagsOut
  estimatesInbreedingVar[[ii]] <- MCMCsummary(coda_samples, round = 4)

  #---------------------------------------------------------------------
  # Computational cost
  #---------------------------------------------------------------------
  par <- rownames(estimatesInbreedingVar[[ii]])
  neef <- estimatesInbreedingVar[[ii]]$n.eff
  rho <- nUseSteps/neef # value 1 is the ideal situation
  cost <- (as.numeric(time_PedInb_mean)*rho)*100
  data_costInbVar[[ii]] <-
    data.frame(par, time_PedInb_mean, time_PedInb_se,
               neef, rho, cost)
  data_costInbVar[[ii]] <- data_costInbVar[[ii]] %>%
    mutate_all(function(x) ifelse(is.infinite(x), NA, x))
  data_costInbVar[[ii]]$Generation <- ii + 4
  ESSInbVar[[ii]] <- data_costInbVar[[ii]] %>%
    summarise_all(
      list(mean = mean, median = median,  sd = sd), na.rm = TRUE
    )
  #---------------------------------------------------------------------
}

data_cost_all_InbVar <- do.call(rbind.data.frame, data_costInbVar)
ESS_all_InbVar <- do.call(rbind.data.frame, ESSInbVar)

save(runJagsModelInbVar, estimatesInbreedingVar, data_cost_all_InbVar,
     ESS_all_InbVar, file = paste0("./JAGS/runJagsOutInb",ii, ".RData"))

#=======================================================================
# Animal Model ignoring inbreeding
#=======================================================================
AnimalModelBUGS

filename <- "model.txt"
writeJags(model = AnimalModelBUGS, filename = filename)

tmp <- runJagsModelNoInb <- estimatesNoInbreeding <-
   data_cost <- ESS <- list()
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
    ped=jags.data[, c(1:3, 6)], noX = TRUE, inbreeding = FALSE,
  )


  # Set Covariates
  tmp[[ii]]$x <- with(jags.data, model.matrix(~year + sex))
  tmp[[ii]]$nB <- ncol(tmp[[ii]]$x)

  ## Set precisions

  tmp[[ii]]$tau2a <- 1 / sigma2a
  tmp[[ii]]$tau2e <- 1 / sigma2e

  #-----------------------------------------------------------------------
  # Jags
  #-----------------------------------------------------------------------
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
  runJagsModelNoInb[[ii]] <- runJagsOut
  estimatesNoInbreeding[[ii]] <- MCMCsummary(coda_samples, round = 4)

  #---------------------------------------------------------------------
  # Computational cost
  #---------------------------------------------------------------------
  par <- rownames(estimatesNoInbreeding[[ii]])
  neef <- estimatesNoInbreeding[[ii]]$n.eff
  rho <- nUseSteps/neef # value 1 is the ideal situation
  cost <- (as.numeric(time_PedInb_mean)*rho)*100
  data_cost[[ii]] <- data.frame(par, time_PedInb_mean, time_PedInb_se,
                                neef, rho, cost)
  data_cost[[ii]] <- data_cost[[ii]] %>%
    mutate_all(function(x) ifelse(is.infinite(x), NA, x))
  data_cost[[ii]]$Generation <- ii + 4
  ESS[[ii]] <- data_cost[[ii]] %>%
    summarise_all(
      list(mean = mean, median = median,  sd = sd), na.rm = TRUE
    )
  #---------------------------------------------------------------------
}

data_cost_all_NoInb <- do.call(rbind.data.frame, data_cost)
ESS_all_NoInb <- do.call(rbind.data.frame, ESS)

save(runJagsModelNoInb, estimatesNoInbreeding, data_cost_all_NoInb, ESS_all_NoInb,
     file = paste0("./JAGS/runJagsOutNoInb",ii, ".RData"))


#-----------------------------------------------------------------------
# Variance unknown
#----------------------------------------------------------------------
AnimalModelBUGSVar

filename <- "model.txt"
writeJags(model = AnimalModelBUGSVar, filename = filename)

tmp <- runJagsModelNoInbVar <- estimatesNoInbreedingVar <-
   data_cost <- ESS <- list()
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
    ped=jags.data[, c(1:3, 6)], noX = TRUE, inbreeding = FALSE,
  )

  # Set Covariates
  tmp[[ii]]$x <- with(jags.data, model.matrix(~year + sex))
  tmp[[ii]]$nB <- ncol(tmp[[ii]]$x)

  #-----------------------------------------------------------------------
  # Jags
  #-----------------------------------------------------------------------
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
    return(list(a = a[[chain]], b = b[[chain]], tau2a = tau2a[[chain]],
                tau2e = tau2e[[chain]],
                .RNG.seed = .RNG.seed, .RNG.name = .RNG.name))
  }
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
  runJagsModelNoInbVar[[ii]] <- runJagsOut
  estimatesNoInbreedingVar[[ii]] <- MCMCsummary(coda_samples, round = 4)

  #---------------------------------------------------------------------
  # Computational cost
  #---------------------------------------------------------------------
  par <- rownames(estimatesNoInbreedingVar[[ii]])
  neef <- estimatesNoInbreedingVar[[ii]]$n.eff
  rho <- nUseSteps/neef # value 1 is the ideal situation
  cost <- (as.numeric(time_PedInb_mean)*rho)*100
  data_cost[[ii]] <- data.frame(par, time_PedInb_mean, time_PedInb_se,
                                neef, rho, cost)
  data_cost[[ii]] <- data_cost[[ii]] %>%
    mutate_all(function(x) ifelse(is.infinite(x), NA, x))
  data_cost[[ii]]$Generation <- ii + 4
  ESS[[ii]] <- data_cost[[ii]] %>%
    summarise_all(
      list(mean = mean, median = median,  sd = sd), na.rm = TRUE
    )
  #---------------------------------------------------------------------
}

data_cost_all_NoInbVar <- do.call(rbind.data.frame, data_cost)
ESS_all_NoInbVar <- do.call(rbind.data.frame, ESS)

save(runJagsModelNoInbVar, estimatesNoInbreedingVar,
     data_cost_all_NoInbVar, ESS_all_NoInbVar,
     file = paste0("./JAGS/runJagsOutNoInbVar",ii, ".RData"))

#-----------------------------------------------------------------------
# Camparing with and without inbreeding
#-----------------------------------------------------------------------
data_comp <- data.frame(NoInb = estimatesNoInbreeding[c(an), 1],
                        Inb = estimatesInbreeding[c(an), 1],
                        Animal = gl(n = length(an),  k = 1,
                                    labels = paste0("Animal ", an)))

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
  Parameter = rownames(estimatesInbreeding)[c(an, bn)],
  Difference = estimatesNoInbreeding[c(an, bn), 1] -
    estimatesInbreeding[c(an, bn), 1],
  Ratio = estimatesNoInbreeding[c(an, bn), 1]/
    estimatesInbreeding[c(an, bn), 1])

data_diff %>%
  kbl(
    caption = "Difference and ratio between parameter estimates (sample of 'an' animals)",
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
AnimalModelBUGSChol

filename <- "model.txt"
writeJags(model = AnimalModelBUGSChol, filename = "model.txt")

tmp <- runJagsModel <- estimatesCholInbreeding <-
   data_cost <- ESS <- list()
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
    ped=jags.data[, c(1:3, 6)], noX = TRUE, Chol = TRUE
  )


  # Set Covariates
  tmp[[ii]]$x <- with(jags.data, model.matrix(~year + sex))
  tmp[[ii]]$nB <- ncol(tmp[[ii]]$x)

  ## Set precisions

  tmp[[ii]]$tau2a <- 1 / sigma2a
  tmp[[ii]]$tau2e <- 1 / sigma2e

#-----------------------------------------------------------------------
# model
#-----------------------------------------------------------------------
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
    time[i] <- difftime(end_time, start_time, units='mins')
    cat(paste("\n", "Generation", ii, "fit number", i, "\n"))
  }
  (time_PedInb_mean <- mean(time))
  (time_PedInb_se <- sd(time)/sqrt(nt))

  # coda samples - MCMC
  coda_samples <- as.mcmc.list(runJagsOut)

  # parameter estimates
  runJagsModel[[ii]] <- runJagsOut
  estimatesCholInbreeding[[ii]] <- MCMCsummary(coda_samples, round = 4)

  #---------------------------------------------------------------------
  # Computational cost
  #---------------------------------------------------------------------
  par <- rownames(estimatesCholInbreeding[[ii]])
  neef <- estimatesCholInbreeding[[ii]]$n.eff
  rho <- nUseSteps/neef # value 1 is the ideal situation
  cost <- (as.numeric(time_PedInb_mean)*rho)*100
  data_cost[[ii]] <- data.frame(par, time_PedInb_mean, time_PedInb_se,
                                neef, rho, cost)
  data_cost[[ii]] <- data_cost[[ii]] %>%
    mutate_all(function(x) ifelse(is.infinite(x), NA, x))
  data_cost[[ii]]$Generation <- ii + 4
  ESS[[ii]] <- data_cost[[ii]] %>%
    summarise_all(
      list(mean = mean, median = median,  sd = sd), na.rm = TRUE
    )
  #---------------------------------------------------------------------
}

data_cost_all <- do.call(rbind.data.frame, data_cost)
ESS_all <- do.call(rbind.data.frame, ESS)

save(runJagsModel, estimatesCholInbreeding, data_cost_all, ESS_all,
     file = paste0("./JAGS/runJagsOutChol",ii, ".RData"))
#=======================================================================
# Chol with variance unknown
#=======================================================================
AnimalModelBUGSCholVar

filename <- "model.txt"
writeJags(model = AnimalModelBUGSCholVar, filename = "model.txt")

tmp <- runJagsModel <- estimatesCholInbreedingVar <-
   data_cost <- ESS <- list()
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
    ped=jags.data[, c(1:3, 6)], noX = TRUE, Chol = TRUE
  )

  # Set Covariates
  tmp[[ii]]$x <- with(jags.data, model.matrix(~year + sex))
  tmp[[ii]]$nB <- ncol(tmp[[ii]]$x)
#-----------------------------------------------------------------------
# model
#-----------------------------------------------------------------------
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
    time[i] <- difftime(end_time, start_time, units='mins')
    cat(paste("\n", "Generation", ii, "fit number", i, "\n"))
  }
  (time_PedInb_mean <- mean(time))
  (time_PedInb_se <- sd(time)/sqrt(nt))

  # coda samples - MCMC
  coda_samples <- as.mcmc.list(runJagsOut)

  # parameter estimates
  runJagsModel[[ii]] <- runJagsOut
  estimatesCholInbreedingVar[[ii]] <- MCMCsummary(coda_samples, round = 4)

  #---------------------------------------------------------------------
  # Computational cost
  #---------------------------------------------------------------------
  par <- rownames(estimatesCholInbreedingVar[[ii]])
  neef <- estimatesCholInbreedingVar[[ii]]$n.eff
  rho <- nUseSteps/neef # value 1 is the ideal situation
  cost <- (as.numeric(time_PedInb_mean)*rho)*100
  data_cost[[ii]] <- data.frame(par, time_PedInb_mean, time_PedInb_se,
                                neef, rho, cost)
  data_cost[[ii]] <- data_cost[[ii]] %>%
    mutate_all(function(x) ifelse(is.infinite(x), NA, x))
  data_cost[[ii]]$Generation <- ii + 4
  ESS[[ii]] <- data_cost[[ii]] %>%
    summarise_all(
      list(mean = mean, median = median,  sd = sd), na.rm = TRUE
    )
  #---------------------------------------------------------------------
}

data_cost_all <- do.call(rbind.data.frame, data_cost)
ESS_all <- do.call(rbind.data.frame, ESS)

save(runJagsModel, estimatesCholInbreedingVar, data_cost_all, ESS_all,
     file = paste0("./JAGS/runJagsOutChol",ii, ".RData"))

#=======================================================================
#-----------------------------------------------------------------------
# Camparing models considering inbreeding -
# model 1: estimates using general inverse matrix
# model 2: estimates using Cholesky decomposition
#-----------------------------------------------------------------------
data_comp <- data.frame(CholInb = estimatesCholInbreeding[c(an), 1],
                        Inb = estimatesInbreeding[c(an), 1],
                        Animal = gl(n = length(an),  k = 1,
                                    labels = paste0("Animal ", an)))

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
  Parameter = rownames(estimatesInbreeding)[c(an, bn)],
  Difference = estimatesCholInbreeding[c(an, bn), 1] -
    estimatesInbreeding[c(an, bn), 1],
  Ratio = estimatesCholInbreeding[c(an, bn), 1]/
    estimatesInbreeding[c(an, bn), 1])

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
set.seed(123)
ord <- sample(1:nrow(example), replace = FALSE)
tmp <- AnimalModelBUGSData(ped=example[ord, ], Chol=TRUE,
                            reorder=TRUE)

## Set precisions
tmp$tau2a <- 1 / sigma2a
tmp$tau2e <- 1 / sigma2e
#-----------------------------------------------------------------------
# Model
#-----------------------------------------------------------------------
AnimalModelBUGSChol
jags_params <- c("a", "b")

time <- NA
for (i in 1:nt){
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
  time[i] <- difftime(end_time, start_time, units='secs')
}

(time_PedOrd_mean <- mean(time))
(time_PedOrd_se <- sd(time)/sqrt(nt))

#save(runJagsOut, file = "runJagsOutOrd.RData")
#save(time_PedOrd_mean, file = "time_PedOrd_mean.RData")
#save(time_PedOrd_se, file = "time_PedOrd_se.RData")

# coda samples - MCMC
coda_samples <- as.mcmc.list(runJagsOut)

# parameter estimates
estimatesOrdInbreedingOrd <- MCMCsummary(coda_samples, round = 4)
estimatesOrdInbreedingOrd

#save(estimatesOrdInbreeding, file = "estimatesOrdInbreeding.RData")
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

## Prepair data
(tmp <- AnimalModelBUGSData(ped=example, ZChol=TRUE))

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
AnimalModelBUGSZStar
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
time_Zstar <- difftime(end_time, start_time, units='secs')

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
time_ZStarSVD <- difftime(end_time, start_time, units='secs')

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
