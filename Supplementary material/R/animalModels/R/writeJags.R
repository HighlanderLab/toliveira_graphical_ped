#=======================================================================
# Writing a JAGS model                                                 #
# Contributors: Thiago P. Oliveira, Ivan Pocrnic,  and Gregor Gorjanc  #
# Written by Thiago P. Oliveira                                        #
                                                                       #
# Functions:                                                           #
                                                                       #
# License: GNU General Public License version 2 (June, 1991) or later  #
# Last Update: 04 Nov 2020                                             #
#=======================================================================
##' @title BUGS to JAGS
##'
##' @description Writes a BUGS model file in a JAGS format
##'
##' @usage
##' writeJags(model, filename, dir, overwrite, digits)
##'
##' @param model JAGS model format to be write onto the hard drive.
##'
##' @param filename character,  the name of the file to write/overwrite.
##'
##' @param dir directory where to write the file. Defaults to uses the
##'   current working directory through \code{getwd}
##'
##' @param overwrite logical, if \code{TRUE} the \code{filename} with
##'   the same name will be overwrited.
##'
##' @param digits Number of digits used in the output.
##'
##' @author Thiago de Paula Oliveira, \email{thiago.oliveira@ed.ac.uk},
##'   Ivan Pocrnic, \email{ivan.pocrnic@roslin.ed.ac.uk}, Gregor
##'   Gorjanc, \email{gregor.gorjanc@gmail.com}
##'
##' @seealso \code{\link{AnimalModelBUGSData}}
##'
##' @keywords JAGS BUGS
##'
##' @examples
##' \dontrun{
##'  AnimalModelBUGSVar <- function()
##'   {
##'    ## Priors - precision parameters
##'    tau2e ~ dgamma(0.001, 0.001)
##'    tau2a ~ dgamma(0.001, 0.001)
##'    sigma2e <- pow(tau2e, -1)
##'    sigma2a <- pow(tau2a, -1)
##'    ## Additive genetic values
##'    for(k in 1:nI) {
##'     a[id[k]] ~ dnorm(pa[id[k]], Xtau2a[id[k]])
##'     pa[id[k]] <- 0.5 * (a[fid[k]] + a[mid[k]])
##'     Xtau2a[id[k]] <- winv[id[k]] * tau2a
##'    }
##'    a[nU] <- 0 # NULL (zero) holder
##'    ## Fixed Effects
##'    for(j in 1:nB) { b[j] ~ dnorm(0, 1.0E-6) }
##'     ## Phenotypes
##'     for(i in 1:nY) {
##'      y[i] ~ dnorm(mu[i], tau2e)
##'      mu[i] <- inprod(x[i, ], b[]) + a[idy[i]]
##'     }
##'   }
##'
##' filename <- "model.txt"
##'
##' writeJags(model = AnimalModelBUGSVar, filename = filename)
##' }
##'
##' @export

writeJags <- function(model,  filename,  dir = getwd(),
                      overwrite =  TRUE,  digits = 5) {
  #
  #---------------------------------------------------------------------
  # General arguments
  #---------------------------------------------------------------------
  file_ext <- "txt"
  get_body <- body(model)
  #---------------------------------------------------------------------
  # getting the model
  #---------------------------------------------------------------------
  model.text <- c("model",
                  replaceScientificNotationR(get_body, digits = digits))
  #---------------------------------------------------------------------
  # Excluding file extension
  #---------------------------------------------------------------------
  sv <- unlist(strsplit(filename, "\\."))
  # Removing possible extension
  if (sv[length(sv)] == file_ext) {
    sv <- paste(sv[-length(sv)],  collapse = ".")
  }
  endFile <- paste0(sv, ".", file_ext)
  #---------------------------------------------------------------------
  # Path
  #---------------------------------------------------------------------
  my_path <- if (is.null(dir)) {
    getwd()
  }else {
    as.character(dir)
  }
  endFile <- file.path(my_path, endFile)
  #---------------------------------------------------------------------
  # Saving
  #---------------------------------------------------------------------
  if (overwrite == TRUE) {
    writeLines(model.text, endFile)
  }else {
    i <- 0
    while (file.exists(endFile) == TRUE) {
      i <- i + 1
      endFile <- paste0(sv, "(", i, ")", ".", file_ext)
    }
    writeLines(model.text, endFile)
  }
}

#=======================================================================
# R2WinBUGS function
# Originally written by Andrew Gelman
#=======================================================================
##' @title Internal Function to write JAGS model
##'
##' @description This is an internally called function used to write
##'   theJAGS model. Modified from Andrew Gelman function.
##'
##' @usage NULL
##'
##' @keywords internal
replaceScientificNotationR <- function(bmodel, digits = 5){
  env <- new.env()
  assign("rSNRidCounter", 0, envir=env)
  replaceID <- function(bmodel, env, digits = 5){
    for(i in seq_along(bmodel)){
      if(length(bmodel[[i]]) == 1){
        if(as.character(bmodel[[i]]) %in% c(":", "[", "[[")) {
          return(bmodel)
        }
        if((typeof(bmodel[[i]]) %in% c("double", "integer")) &&
             ((abs(bmodel[[i]]) < 1e-3) || (abs(bmodel[[i]]) > 1e+4))){
          counter <- get("rSNRidCounter", envir=env) + 1
          assign("rSNRidCounter", counter, envir=env)
          id <- paste("rSNRid", counter, sep="")
          assign(id, formatC(bmodel[[i]], digits=digits, format="E"),
                 envir=env)
          bmodel[[i]] <- id
        }
      } else {
        bmodel[[i]] <- replaceID(bmodel[[i]], env, digits = digits)
      }
    }
    bmodel
  }
  bmodel <- deparse(replaceID(bmodel, env, digits = digits),
                    control = NULL)
  for(i in ls(env)){
    bmodel <- gsub(paste('"', i, '"', sep=''), get(i, envir=env),
                   bmodel, fixed=TRUE)
  }
  bmodel
}
