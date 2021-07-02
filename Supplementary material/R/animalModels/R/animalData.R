#######################################################################
#                                                                     #
# Package: animalData                                                 #
#                                                                     #
# File: animalData                                                    #
# Contains: AnimalModelBUGSData function                              #
#                                                                     #
# Written by Thiago de Paula Oliveira,  Ivan Pocrnic,  Gregor Gorjanc #
# copyright (c) 2021, Oliveira et al.                                 #
#                                                                     #
# First version: 27/01/2021                                           #
# License: GNU General Public License version 2 (June, 1991) or later #
#                                                                     #
#######################################################################

##' @title Prepair data for animal model in BUGS
##'
##' @description Prepair data for pedigree based mixed effects model
##'   (animal model) in BUGS. This function does not handle
##'   separate genetic groups for individual and parental additive
##'   genetic effects! To combine RAM and/or parental and/or multiple
##'   trait models, a user must provide the BUGS models. Only one of
##'   Chol, ZChol, or ZSVD can be used. The later two can not be used in
##'   combination with inbreeding, groups, RAM, and parental options.
##'
##' @usage
##' AnimalModelBUGSData(ped, phen, ncolY, reorder, inbreeding,
##'   inbreedingCoef, groups, RAM, maternal, paternal, perm, uncertain,
##'   pedigreemm, noX, Chol, ZChol, ZSVD, intercept)
##'
##' @param ped an object of class \code{data.frame} or \code{matrix}. A
##'   classic pedigree with three columns with individual, father, and
##'   mother id in numeric format; it can also hold phenotype data after
##'   the third column (see phen), but this is handy only with
##'   non-repeated records. It is assumed that the pedigree is prepaired
##'   so that: \describe{\item{*id:}{range from 1 to n}
##'   \item{\code{NA}:}{is used for an unknown parent} \item{fathers
##'   and mothers:}{appear as individuals (they have their own row) so
##'   that the number of individuals (\code{nI}) is equal to the number
##'   of rows} \item{*parent id:}{are greater than individual id}
##'   \item{*parents:}{we assumed they precede their children otherwise
##'   use \code{reoder=TRUE} (this last assumption is not tested!)}}
##'   Assumptions marked with * are not needed for BUGS, but for R
##'   programs that are used to work with the pedigree in this function.
##'
##' @param phen an object of class \code{data.frame}. Data related to
##'   phenotypes - the first column must be individual id (in numeric
##'   format), then \code{ncolY} columns of phenotypic values and
##'   finally group predictor columns (in numeric format); when
##'   \code{phen} is missing it is assumed that the \code{ped} contains
##'   the corresponding data; currently only group predictors are
##'   supported, but this could be extended quite easily!  Default is
##'   \code{NULL}.
##'
##' @param ncolY numeric, number of columns with phenotypic
##'   values. Default is \code{1}.
##'
##' @param reorder logical, setting \code{reorder=TRUE} causes
##'   reordering of the pedigree (via
##'   \code{\link[MasterBayes]{orderPed}}) for the internal
##'   calculations, while the data is returned in the original order; if
##'   pedigree is not ordered, inbreeding can not be
##'   \code{TRUE}. Default is \code{FALSE}.
##'
##' @param inbreeding logical, account for inbreeding when computing
##'   within family segregation variance coefficients. Default is
##'   \code{TRUE}.
##'
##' @param inbreedingCoef logical, output individual inbreeding
##'   coefficients for phenotyped individuals - not yet
##'   implemented. Default is \code{FALSE}.
##'
##' @param groups logical, does the pedigree contain the so called
##'   genetic groups for unknown parents; when \code{TRUE}, pedigree
##'   should be prepaired as usuall with the following additions (bellow
##'   is an example): \describe{\item{Ex. 1:}{unknown parents are
##'   replaced with phantom parents} \item{Ex. 2:}{genetic groups so
##'   that there are no unknown parents, e.g., no NA, left in the
##'   pedigree} \item{Ex. 3:}{genetic groups must have identifications
##'   from (\code{nI + 1}) to (\code{nI + nG}), where \code{nI} is the
##'   number of individuals in the pedigree (= number of rows) and
##'   \code{}{nG} is the number of genetic groups}}. Default is
##'   \code{FALSE}.
##'
##' @param RAM logical, prepair data for the Reduced Animal Model
##'   (RAM). Default is \code{FALSE}.
##'
##' @param maternal logical, prepair data for maternal model, i.e.,
##'   adding mother id to phenotypic data. Default is \code{FALSE}.
##'
##' @param paternal logical, prepair data for paternal model, i.e.,
##'   adding father id to phenotypic data. Default is \code{FALSE}.
##'
##' @param perm logical, prepair data for individual permanent effect
##'   (this needs to be treated separately, due to BUGS model language -
##'   we can not loop over, say, p[1], p[3], p[10], ...). Default is
##'   \code{FALSE}.
##'
##' @param uncertain logical, there are uncertain parentages in the
##'   pedigree - not yet implemented.  Default is \code{FALSE}.
##'
##' @param pedigreemm logical, use the package \code{\link{pedigreemm}}
##'   or \code{\link{MCMCglmm}} for the calculation of within family
##'   segregation variance coefficients - \code{\link[pedigreemm]{Dmat}}
##'   is faster than \code{\link[MCMCglmm]{inverseA}}, due to less
##'   overhead - \code{\link{pedigreemm}} is still needed for other
##'   options, though.  Default is \code{TRUE}.
##'
##' @param noX logical, there is no "fixed" effects. Default is
##'   \code{FALSE}.
##'
##' @param Chol logical, prepair data for the reparametrized model
##'   (Cholesky in prior) as in Damgaard (2007). Default is
##'   \code{FALSE}.
##'
##' @param ZChol logical, prepair reparametrized design matrix with the
##'   Cholesky decomposition of relationship matrix. Default is
##'   \code{FALSE}.
##'
##' @param ZSVD logical, prepair reparametrized design matrix with the
##'   singular value decomposition of relationship matrix. Default is
##'   \code{FALSE}.
##'
##' @param intercept logical, fit a model with an intercept. Default
##'   is \code{TRUE}.
##'
##' @return return a list of objects.
##'
##' @author Thiago de Paula Oliveira, \email{thiago.oliveira@ed.ac.uk},
##'   Ivan Pocrnic, \email{ivan.pocrnic@roslin.ed.ac.uk}, Gregor Gorjanc,
##'   \email{gregor.gorjanc@gmail.com}
##'
##' @importFrom pedigreemm Dmat relfactor
##'
##' @importFrom MCMCglmm inverseA
##'
##' @importFrom MasterBayes orderPed
##'
##' @importFrom stats complete.cases model.matrix
##'
##' @importFrom utils packageDescription
##'
##' @importFrom Matrix Diagonal t crossprod
##'
##' @seealso \code{\link{writeJags}}
##'
##' @references Damgaard, L. H. Technical note: How to use Winbugs to
##'   draw inferences in animal models. \emph{Journal of Animal
##'   Science}, v. 85, n. 6, p. 1363-1368, 2007
##'
##' @keywords BUGS  Mixed-Effects-Model  Pedigree
##'
##' @examples
##' ## Preparing the data
##' example <- data.frame(
##'  individual = c(  1,   2,   3,   4,   5,   6,   7,   8,   9),
##'  father = c( NA,  NA,   2,   2,   4,   2,   5,  NA,   7),
##'  mother = c( NA,  NA,   1,  NA,   3,   3,   6,  NA,   8),
##'  phenotype = c( NA, 105,  98, 101, 106,  93,  NA,  NA, 109),
##'  group = factor(c( NA,   1,   1,   2,   2,   2,  NA,  NA,   1))
##' )
##' tmp <- AnimalModelBUGSData(ped=example)
##' print(tmp)
##'
##' @export

AnimalModelBUGSData <- function(ped, phen=NULL, ncolY=1, reorder=FALSE,
                                inbreeding=TRUE, inbreedingCoef=FALSE,
                                groups=FALSE, RAM=FALSE, maternal=FALSE,
                                paternal=FALSE, perm=FALSE,
                                uncertain=FALSE, pedigreemm=TRUE,
                                noX=FALSE, Chol=FALSE, ZChol=FALSE,
                                ZSVD=FALSE,  intercept = TRUE)
{
  ## --- Some tests ---

  ped <- as.data.frame(ped)

  if(inbreedingCoef) stop("'inbreedingCoef' not yet implemented")
  if(uncertain) stop("'uncertain' not yet implemented")

  if(is.null(phen)) phen <- ped[, c(1, 4:ncol(ped))]

  test <- !all(sapply(ped[, 1:3], is.numeric))
  if(test) stop("identifications in 'ped' must be in numeric format")

  test <- !is.numeric(phen[, 1])
  if(test) stop("identifications in 'phen' must be in numeric format")

  test <- (min(ped[, 1]) < 1) || (max(ped[, 1]) < nrow(ped))
  if(test) stop("indentifications must range between 1 and n")

  if(!groups) {
    tmp <- c(ped[, 2], ped[, 3])
    test <- any(!(tmp[!is.na(tmp)] %in% ped[, 1]))
    if(test) stop("known fathers and mothers must appear as individuals in 'ped'")

    no.father <- is.na(ped[, 2])
    no.mother <- is.na(ped[, 3])
    if(!reorder &
      any(ped[!no.father, 2] >= ped[!no.father, 1]) ||
      any(ped[!no.mother, 3] >= ped[!no.mother, 1]))
        stop("father and mother id must be greater than individual id")
  }

  test <- any(!(phen[, 1] %in% ped[, 1]))
  if(test) stop("individuals in 'phen' must have an entry in 'ped'")

  test <- (Chol + ZChol + ZSVD) > 1
  if(test) stop("only one reparametrization ('Chol', 'ZChol', or 'ZSVD') can be used at a time")

  test <- (ZChol | ZSVD) & any(c(!inbreeding, groups, RAM, maternal, paternal))
  if(test) stop("'ZChol' and 'ZSVD' can not be used with 'inbreeding', 'groups', 'RAM', and parental options")

  ## --- Phenotype data ---

  ## Columns:
  ## 1                - id
  ## 2:(1 + ncolY)    - phenotypic values
  ## (2 + ncolY):ncol - predictors

  if(!noX) { # Model will have some group predictors and/or covariates
    ## Keep the rows that have all predictors known
    ##  --> missing phenotypes can be automatically imputed by BUGS
    phen <- phen[complete.cases(phen[, c(1, (1 + ncolY + 1):ncol(phen))]), ]
  } else {   # Model will not have group predictors and/or covariates
    ## Keep all the rows that have at least one phenotype known
    ##  --> missing phenotypes can be automatically imputed by BUGS
    tmp <- !is.na(phen[, 1:(1 + ncolY)])
    if(ncolY > 1) tmp <- as.logical(rowSums(tmp))
    phen <- phen[tmp[, 1], ]
  }
  nY <- nrow(phen)

  ## --- Pedigree related stuff ---

  ## Reorder
  if(reorder) {
    orderOrig <- as.character(ped[, 1])
    ped <- MasterBayes::orderPed(ped=ped)
    ## ToDo: what about pedigree::orderPed
  }

  if(!ZChol & !ZSVD) {

    nI <- nrow(ped)
    ## Within family segregation variance coefficients

    if(groups) {  ## Otherwise inverseA() complains that some parents do
      ped2 <- ped ## not appear as individuals
      ped[ped[, 2] > nI, 2] <- NA
      ped[ped[, 3] > nI, 3] <- NA
    }
    no.father <- is.na(ped[, 2])
    no.mother <- is.na(ped[, 3])

    if(inbreeding) {
      ## Due to buglet in pedigreemm 0.2.3
      ver <- package_version(packageDescription("pedigreemm")$Version) >= "0.2.4"
      test <- (reorder + ver) %in% c(0, 2)
      if(pedigreemm & test) { ## Using Dmat() from pedigreemm
        w <- pedigreemm::Dmat(ped=pedigreemm::pedigree(sire=ped[, 2], dam=ped[, 3], label=ped[, 1]))
      } else {                ## Using inverseA() from MCMCglmm
        w <- MCMCglmm::inverseA(pedigree=ped[, 1:3])$dii
      }
    } else {
      w <- rep(0, times=nI)
      w[no.father & no.mother] <- 1
      w[(no.father & !no.mother) | (!no.father & no.mother)] <- 0.75
      w[w < 0.75] <- 0.5
    }

    if(groups) {
      ## Add 1/4m, where m is number of base parents; see Quass (1988) page 1340 bottom left
      ## However, original w is computed with possibility of missing parents, but parents are
      ## with groups assumed known!
      ## With groups A = inv(I-P) (1/4m + var(w)) inv(I-t(P))
      ## No parents known w = 1, but with groups parents are known so set w = 1 - 1/4 - 1/4 = 1/2
      w[no.father & no.mother] <- 0.5

      ## One parent known w = 1 - 1/4(1 + F[parent]), but with groups both parents are known so
      ## set w = 1 - 1/4(1 + F[parent] - 1/4
      tmp <- no.father & !no.mother
      w[tmp] <- w[tmp] - 1/4
      tmp <- !no.father &  no.mother
      w[tmp] <- w[tmp] - 1/4

      ## Now add 1/4m
      w <- w + 0.25 * (no.father + no.mother)
      ped <- ped2
      no.father <- is.na(ped[, 2])
      no.mother <- is.na(ped[, 3])
    }

    winv <- 1 / w
    names(winv) <- ped[, 1]
    ## Making sure that winv is sorted by id as 1:nI - we need this for RAM,
    ## so it is done for all models for completeness
    winv <- winv[as.character(1:nI)]
    names(winv) <- NULL

    ## Unknown parents / genetic groups
    nU <- nI + 1
    if(!groups) {
      ## Modify unknown parents to nI + 1 (this will be a NULL (zero) holder in BUGS)
      ped[no.father, 2] <- nU
      ped[no.mother, 3] <- nU
    } else {
      ## Test
      if(any(no.father) | any(no.mother)) stop("When 'groups=TRUE', all unknown parents must be assigned to genetic groups!")
      ## Number of genetic groups
      nG <- max(ped[, 2:3]) ## no need for na.rm=TRUE in max(), because it is already tested above
      if(nG == nU) warning("It appears that there is only one one genetic group!")
    }

  } else {

    ## Create design matrix Z for reparametrization - pedigree will be infiltrated into this matrix
    tmp <- factor(phen[, 1]) ## only those individuals that have records!
    nI <- nlevels(tmp)
    Z <- model.matrix(~ tmp - 1) #changed 21 Jan 2021 - as.matrix
    U <- relfactor(ped=pedigreemm::pedigree(sire=ped[, 2], dam=ped[, 3], label=ped[, 1]), labs=levels(tmp))

    if(ZChol) { ## Parametrization via the Cholesky decomposition of A (Quaas, 1984)
      ZStar <- as.matrix(Z %*% Matrix::t(U))
      ## Remove some smallish values
      ZStar[ZStar < 0.00000001] <- 0
    }

    if(ZSVD) {  ## Parametrization via the singular value decomposition of A (Waldmann et al., 2008)
      A <- Matrix::crossprod(U)            ## first bottle-neck (the A matrix)
      svdUD <- svd(A)[c("u", "d")] ## second bottle-neck (the computations and the svdUD object)
      ZStar <- as.matrix(Z %*% svdUD$u %*% Diagonal(x=sqrt(svdUD$d)))
    }

  }

  ## Reorder pedigree back to the original order
  if(reorder) {
    rownames(ped) <- ped[, 1]
    ped <- ped[orderOrig, ]
    rownames(ped) <- NULL
  }

  ## --- Modifications for the Reduced Animal Model (RAM) ---

  if(RAM) {
    ## Add non-parent indicator in the data
    phen <- cbind(phen, !is.parent(ped=ped, id=phen[, 1]))
    ## Add non-parent indicator in the pedigree
    ped <- cbind(ped, !is.parent(ped=ped))
    ## Remove non-parents from the pedigree
    pedFull <- ped
    ped <- ped[!ped[, ncol(ped)], ] ## ncol(ped) is parent indicator
    nI <- nrow(ped)
    ## nU or nG stays as before, because identifications can be larger than
    ## nI, i.e., it is not necessary that non-parents have the last (maximal)
    ## id numbers; BUGS seems to work without problems if there are "gaps" in
    ## the parameter vectors, e.g., a[1], a[2], a[4] works fine
  }

  ## --- Return ---

  ## We assume here the following order of columns
  ##
  ## Pedigree:
  ##  ped[, 1] - individual id
  ##  ped[, 2] - father id
  ##  ped[, 3] - mother id
  ##
  ## Data:
  ##  phen[, 1]                      - individual id
  ##  phen[, 2:(1 + ncolY)]          - phenotypic variable(s)
  ##  phen[, (2 + ncolY):(ncol-RAM)] - group/predictor variable
  ##  if(RAM) phen[, ncol]           - non-parent indicator

  ret <- list(nI=nI, nY=nY, y=as.matrix(phen[, 2:(1 + ncolY)])[, , drop=TRUE])

  if(!ZChol & !ZSVD) {
    ret$nU <- nU
    ret$id=ped[, 1]
    ret$fid=ped[, 2]
    ret$mid=ped[, 3]
    ret$winv=winv
    ret$idy=phen[, 1]
  }

  if(ZChol | ZSVD) {
    ret$ZStar <- ZStar
    ret$id <- phen[, 1]
    if(ZChol) ret$fact <- U
    if(ZSVD) ret$fact <- svdUD
  }

  if(groups) ret$nG <- nG

  if(RAM) {
    ret$fidy <- father(ped=pedFull, id=phen[, 1])
    ret$midy <- mother(ped=pedFull, id=phen[, 1])
    ret$npI <- as.integer(phen[, ncol(phen)])
    phen[, ncol(phen)] <- NULL ## no need to care about this column (non-parent indicator) anymore
    ret$npI2 <- apply(X=cbind(ret$idy, ifelse(groups, ret$nG + 1, ret$nU) * ret$npI), MARGIN=1, FUN=max)
  }

  if(maternal) {
    ret$midy <- mother(ped=ped, id=phen[, 1])
  }

  if(paternal) {
    ret$fidy <- father(ped=ped, id=phen[, 1])
  }

  if(perm) {
    ret$idp <- unique(phen[, 1])
    ret$nP <- length(ret$idp)
  }

  if(ncolY > 1) {
    ret$nT <- ncolY
    colnames(ret$y) <- NULL
  }

  if(Chol) {
    ret$wsqr <- sqrt(1 / ret$winv)
    ret$winv <- NULL
  }

  if(!noX) {
    data_phen <- data.frame(phen[, (2 + ncolY):ncol(phen)])
    colnames(data_phen) <- colnames(phen)[(2 + ncolY):ncol(phen)]
    if (intercept == FALSE) {
      ret$x=model.matrix(~ . -1, data = data_phen)
    }else {
      ret$x=model.matrix(~ ., data = data_phen)
    }
    ret$nB=ncol(ret$x)
    names(ret$nB) <- NULL
  }

  ret

  ## Pedigree related:
  ##  nI - number of individuals in the pedigree
  ##  nU - upper bound index/pointer for individuals
  ##  id - individual id
  ##  fid - father id for each individual
  ##  mid - mother ...
  ##  winv - inverse of within family segregation variance coefficients
  ##         (sorted by id as 1:nI)
  ##  wsqr - (for Chol) square root of within family segregation variance
  ##         coefficients (sorted by id as 1:nI)
  ##  nG - (nG - nU + 1) = number of genetic groups
  ##
  ## Data related:
  ##  nY - number of phenotypic records
  ##  idy - individual id
  ##  fidy - (for RAM or paternal) father id for each individual
  ##  midy - (for RAM or maternal) mother ...
  ##  y - phenotypes
  ##  nB - number of location parameters (possibly many)
  ##  x - group/predictor variable(s)
  ##  F - inbreeding coefficient that can be used as a predictor in the model
  ##  npI (for RAM) - non-parent indicator (1 - non-parent, 0 - parent)
  ##  npI2 (for RAM) - parent id or pointer to the NULL (zero) holder in
  ##                   BUGS, i.e., nU (basic RAM) or nG + 1 (RAM with GG)
  ##  ZStar (for ZChol or ZSVD) - reparametrized design matrix for additive
  ##                              genetic values (Z)
  ##  fact  (for ZChol or ZSVD) - factorized relationship matrix
  ##                               - (ZChol) right Cholesky factor (U)
  ##                               - (ZSVD) matrices U and D from svd()
  ##  idp - individuals with phenotypes - for permanent effect
  ##  nP - number of permanent effect levels
}

#=======================================================================
# General methdos
#=======================================================================
##' @title Internal Function to check parents
##'
##' @description This is an internally called function used to check
##'   parents.
##'
##' @usage NULL
##'
##' @keywords internal

is.parent <- function(ped, id)
{
  ## Utility function that returns TRUE or FALSE if individual appears as a
  ## parent or not, respectively
  ##
  ## ped - data.frame, a pedigree holding three columns (individual, father, and mother)
  ## id - individuals for which a test is done, default is ped[, 1]

  if(missing(id)) id <- ped[, 1]

  ## Convert to character to avoid possible problems with combining factors
  id %in% c(as.character(ped[, 2]), as.character(ped[, 3]))
}


##' @title Internal Function to check parents
##'
##' @description This is an internally called function used to check
##'   parents.
##'
##' @usage NULL
##'
##' @keywords internal

father <- function(ped, id)
{
  ## Utility function that returns a father of individual id from the pedigree.
  ##
  ## ped - data.frame, a pedigree holding three columns (individual, father, and mother)
  ## id - individuals for which a father is returned, default is ped[, 1]
  ##
  ## ToDo: add an option to treat unknown parents (=NA) differently, say adding
  ##       phantom parents

  missId <- FALSE
  if(missing(id)) { missId <- TRUE; id <- ped[, 1] }

  ## Test
  if(!missId) {
    test <- any(!(id %in% ped[, 1]))
    if(test) stop("individuals in 'id' must be present in ped[, 1]")
  }
  test <- ped[, 1] %in% id
  ped[test, 2]
}

##' @title Internal Function to check parents
##'
##' @description This is an internally called function used to check
##'   parents.
##'
##' @usage NULL
##'
##' @keywords internal

mother <- function(ped, id)
{
  ## Utility function that returns a mother of individual id from the pedigree.
  ##
  ## ped - data.frame, a pedigree holding three columns (individual, father, and mother)
  ## id - individuals for which a mother is returned, default is ped[, 1]
  ##
  ## ToDo: add an option to treat unknown parents (=NA) differently, say adding
  ##       phantom parents

  missId <- FALSE
  if(missing(id)) { missId <- TRUE; id <- ped[, 1] }

  ## Test
  if(!missId) {
    test <- any(!(id %in% ped[, 1]))
    if(test) stop("individuals in 'id' must be present in ped[, 1]")
  }
  test <- ped[, 1] %in% id
  ped[test, 3]
}

##' @title Internal Function to extend pedigree
##'
##' @usage NULL
##'
##' @keywords internal

ped2pedGG <- function(ped)
{
  ## Extend pedigree so that all unknown parents become known -> phantom parents
  ##
  ## ped - data.frame, having three integer columns: id, fid, and mid all
  ##       in range from 1:n, where n is number of rows; unknown parents
  ##       are denoted with NA
  ##
  ## Result is an augmented and 1:n recoded pedigree

  nI <- nrow(ped)
  test <- !sapply(ped, is.integer)
  if(any(test)) stop("all columns must have integer values!")
  no.father <- is.na(ped$father)
  no.mother <- is.na(ped$mother)
  nF <- sum(no.father)
  nM <- sum(no.mother)
  ped <- ped + (nF + nM)
  ped[no.father, 2] <- 1:nF
  ped[no.mother, 3] <- (1 + nF):(nF + nM)
  ped2 <- data.frame(1:(nF + nM),
                     NA,
                     NA)
  colnames(ped2) <- colnames(ped)
  ped <- rbind(ped2, ped)
  ped
}

##' @title Internal Function to create pedigree uncertain
##'
##' @usage NULL
##'
##' @importFrom gdata NAToUnknown trim
##'
##' @keywords internal

dat2pedigree_uncertain <- function(x, sep="|")
{
  ## Create pedigree_uncertain object
  ##
  ## x - data.frame, having five columns (the rest is ignored): (1) individual,
  ##     (2) father and (3) mother identification and (4, 5) two columns with
  ##     probabilities of parentage for father(s) and mother(s); all values must
  ##     be characters, though the content must be numeric; in case of uncertain
  ##     parentage, the values must be like "2|3" for identifications and
  ##     "0.67|0.33" for probabilities; unknown parents must be presented with NA;
  ##     for parents with certain parentage the probabilities can be NA
  ##
  ## sep - character, used to limit multiple parents and their probabilities
  ##       of parentage
  ##
  ## This code assumes that identifications are already in 1:2 form and that
  ## parents preceede children!
  ##
  ## ToDo: A LOT!
  ##  - make sure all parents appear as individuals
  ##  - other tests
  ##  - if probs are missing and there is more than one possible parent, say
  ##    father, then assume that probs are 1/n for each of n fathers
  ##  - do we need to make sure that uncertain parents must be ordered, i.e.,
  ##    1, 2, 3 or can in it be like 2, 3, 1 --> how does this behave in the
  ##    functions that use this object?

  ## --- Setup ---

  ## To make strsplit happy :( and making sure there are no space
  x <- trim(s=NAToUnknown(x=x[, 1:5], unknown=""))

  ## Test if any "uncertain" string starts or ends with a separator
  testBegEnd <- function(x, string) {
    test <- grep(pattern=paste("^\\", sep, sep=""), x=x)
    if(length(test) > 1) stop(paste(string, "must not start with the separator"))
    test <- grep(pattern=paste("\\", sep, "$", sep=""), x=x)
    if(length(test) > 1) stop(paste(string, "must not end with the separator"))
  }
  testBegEnd(x=x[, 2], string="father id")
  testBegEnd(x=x[, 3], string="mother id")
  testBegEnd(x=x[, 4], string="father prob")
  testBegEnd(x=x[, 5], string="mother prob")

  ## --- Code ---

  ped <- vector(mode="list", length=nrow(x))
  names(ped) <- x[, 1]
  tmp <- vector(mode="list", length=4)
  ## fid, mid, fp, mp

  splitVal <- function(x, sep) {
    ## Split a value by sep and make sure that we get out numeric
    ## which should be 0 in case of NA input
    x <- as.numeric(unlist(strsplit(x, split=sep, fixed=TRUE)))
    if(length(x) < 1) x <- NA
    x
  }

  ## Fill the structure
  for(i in 1:nrow(x)) {
    ped[[i]] <- tmp
    ped[[i]][] <- lapply(x[i, 2:5], splitVal, sep=sep)
    ## The following code assumes that if NA is the first value, then
    ##   the parentage is known with certainity
    if(is.na(ped[[i]][[3]][1])) ped[[i]][[3]] <- 1
    if(is.na(ped[[i]][[4]][1])) ped[[i]][[4]] <- 1
  }

  ## --- Return ---

  ## a list with n components that are lists each with 4 components (vectors)
  ##  - fid - father id
  ##  - mid - mother id
  ##  - fp  - father probability of parentage
  ##  - mp  - mother probability of parentage

  class(ped) <- c("list", "pedigree_uncertain")
  ped
}


##' @title Internal Function to create pedigree uncertain
##'
##' @usage NULL
##'
##' @importFrom gdata NAToUnknown trim
##'
##' @keywords internal

pedigree_uncertain2BUGS <- function(ped, matrix=TRUE)
{
  ## Extract needed data for BUGS from pedigree_uncertain object
  ##
  ## ped - pedigree_uncertain
  ## matrix - logical, return a matrix or (ragged) array

  ## --- Setup ---

  if(!"pedigree_uncertain" %in% class(ped)) stop("'ped' must be of a class 'pedigree_uncertain'")
  n <- length(ped)

  ## Change NA values to n + 1
  ## Assume that if [1] is NA, the parents are unknown
  ped <- lapply(ped, function(x) lapply(x, function(y) if(is.na(y[1])) { n + 1 } else {y}))

  ## --- Code ---

  ## Number of parents
  nC <- sapply(ped, function(x) length(c(x[[1]], x[[2]])))
  names(nC) <- NULL

  if(matrix) {
    ## pid (parent id) and pc (parentage contributions) matrices
    pid <- pc <- matrix(data=0, nrow=n, ncol=max(nC))

    i <- 1
    while(i <= n) {
      pid[i, 1:nC[i]] <- c(ped[[i]][[1]], ped[[i]][[2]]) ## fid, mid
       pc[i, 1:nC[i]] <- c(ped[[i]][[3]], ped[[i]][[4]]) ## probs
      i <- i + 1
    }

    ## --- Return ---

    ## A list with three components:
    ##  - nC  an n*1 vector with number of parent contributions
    ##  - pid an n*max(nC) matrix with parent id
    ##  - pc  an n*max(nC) matrix with parent contributions (must be divided by 2!)

    ret <- list(nC=nC, pid=pid, pc=pc / 2)

  } else {

    ## pid (parent id) and pc (parentage contributions) matrices
    pid <- pc <- vector(mode="integer", length=sum(nC))
    nC2 <- cumsum(nC)
    nC1 <- nC2 - nC + 1

    i <- 1
    while(i <= n) {
      pid[nC1[i]:nC2[i]] <- c(ped[[i]][[1]], ped[[i]][[2]]) ## fid, mid
       pc[nC1[i]:nC2[i]] <- c(ped[[i]][[3]], ped[[i]][[4]]) ## probs
      i <- i + 1
    }

    ## --- Return ---

    ## A list with four components:
    ##  - nC1 an n*1 vector of indexes that start parent contributions of individual i
    ##  - nC2 an n*1 vector of indexes that end   parent contributions of individual i
    ##  - pid an sum(nC)*1 vector with parent id
    ##  - pc  an sum(nC)*1 vector with parent contributions (must be divided by 2!)

    ret <- list(nC1=nC1, nC2=nC2, pid=pid, pc=pc / 2)
  }

  ret
}


##' @title Internal Function to Compute Average Numerator Relationship
##'   Matrix
##'
##' @usage NULL
##'
##' @importFrom gdata NAToUnknown trim
##'
##' @keywords internal

ANRM <- function(ped)
{
  ## Compute Average Numerator Relationship Matrix
  ##  following Perez-Enciso & Fernando (1992)
  ##
  ## ped - pedigree_uncertain
  ##

  ## --- Setup ---

  if(!"pedigree_uncertain" %in% class(ped)) stop("'ped' must be of a class 'pedigree_uncertain'")
  n <- length(ped)

  ## Change NA values to 0
  ## Assume that if [1] is NA, the parents are unknown
  ped <- lapply(ped, function(x) lapply(x, function(y) if(is.na(y[1])) { 0 } else {y}))

  ## Set A
  A <- matrix(data=0, nrow=n, ncol=n)

  ## --- Code ---

  i <- 1
  while(i <= n) {

    ## Off-diagonals
    j <- 1
    while(i > j) {
      A[i, j] <- A[j, i] <- 0.5 * (sum(ped[[i]][[3]] * A[j, ped[[i]][[1]]]) +
                                   sum(ped[[i]][[4]] * A[j, ped[[i]][[2]]]))
      j <- j + 1
    }

    ## Diagonals
    tmp <- A[ped[[i]][[1]], ped[[i]][[2]], drop=FALSE]
    if(nrow(tmp) < 1 || ncol(tmp) < 1) tmp <- 0
    A[i, i] <- 1 + 0.5 * sum(outer(ped[[i]][[3]], ped[[i]][[4]]) * tmp)

    i <- i + 1
  }

  ## --- Return ---

  ## An n*n matrix with relationship coefficients

  A
}

##' @title Internal Function to Compute within family segregation
##'   Matrix
##'
##' @usage NULL
##'
##' @importFrom gdata NAToUnknown trim
##'
##' @keywords internal

ANRMW <- function(ped, A=NULL)
{
  ## Compute within family segregation (Mendelian sampling) variance coefficients
  ##  following Perez-Enciso & Fernando (1992)
  ##
  ## ped - pedigree_uncertain
  ## A - matrix, result of ANRM(); if NULL ANRM() is called

  ## --- Setup ---

  if(!"pedigree_uncertain" %in% class(ped)) stop("'ped' must be of a class 'pedigree_uncertain'")
  n <- length(ped)

  if(is.null(A)) A <- ANRM(ped=ped)

  ## Change NA values to 0
  ## Assume that if [1] is NA, the parents are unknown
  ped <- lapply(ped, function(x) lapply(x, function(y) if(is.na(y[1])) { 0 } else {y}))

  w <- vector(mode="numeric", length=n)

  ## --- Code ---

  i <- 1
  while(i <= n) {

    tmp1 <- A[ped[[i]][[1]], ped[[i]][[1]], drop=FALSE]
    tmp2 <- A[ped[[i]][[2]], ped[[i]][[2]], drop=FALSE]
    if(nrow(tmp1) < 1) tmp1 <- 0
    if(nrow(tmp2) < 1) tmp2 <- 0
    w[i] <- 1 - 0.25 * sum(outer(ped[[i]][[3]], ped[[i]][[3]]) * tmp1) -
                0.25 * sum(outer(ped[[i]][[4]], ped[[i]][[4]]) * tmp2)

    i <- i + 1
  }

  ## --- Return ---

  ## An n*1 vector with family segregation (Mendelian sampling) variance coefficients

  w
}

##' @title Internal Function to Compute diagonals of L and W Matrix
##'
##' @usage NULL
##'
##' @importFrom gdata NAToUnknown trim
##'
##' @keywords internal

FamulaW <- function(ped)
{
  ## Compute diagonals of L and W, where inv(A) =  following (Famula, 1992)
  ##
  ## ped - pedigree_uncertain
  ##
  ## ToDo:
  ##  - this function does not work as it should - find a bug!

  ## --- Setup ---

  if(!"pedigree_uncertain" %in% class(ped)) stop("'ped' must be of a class 'pedigree_uncertain'")
  n <- length(ped)

  ## Change NA values to 0
  ## Assume that if [1] is NA, the parents are unknown
  ped <- lapply(ped, function(x) lapply(x, function(y) if(is.na(y[1])) { 0 } else {y}))

  ## --- Code ---

  ## Working vectors
  w <- v <- u <- vector(mode="numeric", length=n)

  ## Step 1
  i <- 1
  while(i <= n) {

    ## Step 2
    ## Assume that if [1] is NA, the parents are unknown
    test <- sum(is.na(c(ped[[i]][[1]][1], ped[[i]][[2]][1])))
    if(test < 1) { ## Both parents unknown
      v[i] <- 1
      u[i] <- u[i] + 1
    } else {       ## At least one parent known
      v[i] <- sqrt(1 - u[i] + 0.5 * w[i])
      u[i] <- u[i] + v[i] * v[i]
    }

    ## Step 3
    l <- i + 1
    while(l <= n) {

      tmp <- ped[[l]][[1]] >= i
      tf <- sum(v[ped[[l]][[1]][tmp]] * ped[[l]][[3]][tmp])
      tmp <- ped[[l]][[2]] >= i
      tm <- sum(v[ped[[l]][[2]][tmp]] * ped[[l]][[4]][tmp])
      t <- tf + tm
      u[l] <- u[l] + 0.25 * t * t
      w[l] <- w[l] + tf * tm
      v[l] <- 0.5 * t

      l <- l + 1
    }

    i <- i + 1
  }

  ## --- Return ---

  ## A list with two vectors:
  ##  - u diagonal elements of A matrix
  ##  - v diagonal elements of L matrix, where A=Lt(L)

  list(u=u, v=v)
}
