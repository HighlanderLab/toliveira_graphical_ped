## Functions

RecSysMaleFemale <- function(datafile, popname) {
  datafile = rbind(datafile,
                   tibble(generation = generation,
                          ind        = popname@id,
                          father     = popname@father,
                          mother     = popname@mother,
                          sex        = popname@sex,
                          year       = year,
                          pheno      = popname@pheno,
                          tbv        = popname@gv))
}


meanP.2 = function(pop, na.rm = FALSE){   
  colMeans(pop@pheno, na.rm = na.rm)
}

PullSumm = function(datafile, popname){
  gePa = genParam(popname)
  meanP = ifelse(all(is.na(popname@pheno)), NA, meanP.2(popname, na.rm = TRUE))
  datafile = rbind(datafile,
                   tibble(generation = generation,
                          sex        = popname@sex,
                          year       = year,
                          MeanPheno  = meanP,
                          MeanGeno   = gePa$mu,
                          MeanA      = mean(gePa$gv_a),
                          VarA       = gePa$varA,
                          GenicVA    = gePa$genicVarA))
  return(datafile)
}

