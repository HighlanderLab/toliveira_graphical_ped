% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/animalData.R
\name{AnimalModelBUGSData}
\alias{AnimalModelBUGSData}
\title{Prepair data for animal model in BUGS}
\usage{
AnimalModelBUGSData(ped, phen, ncolY, reorder, inbreeding,
  inbreedingCoef, groups, RAM, maternal, paternal, perm, uncertain,
  pedigreemm, noX, Chol, ZChol, ZSVD, intercept)
}
\arguments{
\item{ped}{an object of class \code{data.frame} or \code{matrix}. A
classic pedigree with three columns with individual, father, and
mother id in numeric format; it can also hold phenotype data after
the third column (see phen), but this is handy only with
non-repeated records. It is assumed that the pedigree is prepaired
so that: \describe{\item{*id:}{range from 1 to n}
\item{\code{NA}:}{is used for an unknown parent} \item{fathers
and mothers:}{appear as individuals (they have their own row) so
that the number of individuals (\code{nI}) is equal to the number
of rows} \item{*parent id:}{are greater than individual id}
\item{*parents:}{we assumed they precede their children otherwise
use \code{reoder=TRUE} (this last assumption is not tested!)}}
Assumptions marked with * are not needed for BUGS, but for R
programs that are used to work with the pedigree in this function.}

\item{phen}{an object of class \code{data.frame}. Data related to
phenotypes - the first column must be individual id (in numeric
format), then \code{ncolY} columns of phenotypic values and
finally group predictor columns (in numeric format); when
\code{phen} is missing it is assumed that the \code{ped} contains
the corresponding data; currently only group predictors are
supported, but this could be extended quite easily!  Default is
\code{NULL}.}

\item{ncolY}{numeric, number of columns with phenotypic
values. Default is \code{1}.}

\item{reorder}{logical, setting \code{reorder=TRUE} causes
reordering of the pedigree (via
\code{\link[MasterBayes]{orderPed}}) for the internal
calculations, while the data is returned in the original order; if
pedigree is not ordered, inbreeding can not be
\code{TRUE}. Default is \code{FALSE}.}

\item{inbreeding}{logical, account for inbreeding when computing
within family segregation variance coefficients. Default is
\code{TRUE}.}

\item{inbreedingCoef}{logical, output individual inbreeding
coefficients for phenotyped individuals - not yet
implemented. Default is \code{FALSE}.}

\item{groups}{logical, does the pedigree contain the so called
genetic groups for unknown parents; when \code{TRUE}, pedigree
should be prepaired as usuall with the following additions (bellow
is an example): \describe{\item{Ex. 1:}{unknown parents are
replaced with phantom parents} \item{Ex. 2:}{genetic groups so
that there are no unknown parents, e.g., no NA, left in the
pedigree} \item{Ex. 3:}{genetic groups must have identifications
from (\code{nI + 1}) to (\code{nI + nG}), where \code{nI} is the
number of individuals in the pedigree (= number of rows) and
\code{}{nG} is the number of genetic groups}}. Default is
\code{FALSE}.}

\item{RAM}{logical, prepair data for the Reduced Animal Model
(RAM). Default is \code{FALSE}.}

\item{maternal}{logical, prepair data for maternal model, i.e.,
adding mother id to phenotypic data. Default is \code{FALSE}.}

\item{paternal}{logical, prepair data for paternal model, i.e.,
adding father id to phenotypic data. Default is \code{FALSE}.}

\item{perm}{logical, prepair data for individual permanent effect
(this needs to be treated separately, due to BUGS model language -
we can not loop over, say, p[1], p[3], p[10], ...). Default is
\code{FALSE}.}

\item{uncertain}{logical, there are uncertain parentages in the
pedigree - not yet implemented.  Default is \code{FALSE}.}

\item{pedigreemm}{logical, use the package \code{\link{pedigreemm}}
or \code{\link{MCMCglmm}} for the calculation of within family
segregation variance coefficients - \code{\link[pedigreemm]{Dmat}}
is faster than \code{\link[MCMCglmm]{inverseA}}, due to less
overhead - \code{\link{pedigreemm}} is still needed for other
options, though.  Default is \code{TRUE}.}

\item{noX}{logical, there is no "fixed" effects. Default is
\code{FALSE}.}

\item{Chol}{logical, prepair data for the reparametrized model
(Cholesky in prior) as in Damgaard (2007). Default is
\code{FALSE}.}

\item{ZChol}{logical, prepair reparametrized design matrix with the
Cholesky decomposition of relationship matrix. Default is
\code{FALSE}.}

\item{ZSVD}{logical, prepair reparametrized design matrix with the
singular value decomposition of relationship matrix. Default is
\code{FALSE}.}

\item{intercept}{logical, fit a model with an intercept. Default
is \code{TRUE}.}
}
\value{
return a list of objects.
}
\description{
Prepair data for pedigree based mixed effects model
  (animal model) in BUGS. This function does not handle
  separate genetic groups for individual and parental additive
  genetic effects! To combine RAM and/or parental and/or multiple
  trait models, a user must provide the BUGS models. Only one of
  Chol, ZChol, or ZSVD can be used. The later two can not be used in
  combination with inbreeding, groups, RAM, and parental options.
}
\examples{
## Preparing the data
example <- data.frame(
 individual = c(  1,   2,   3,   4,   5,   6,   7,   8,   9),
 father = c( NA,  NA,   2,   2,   4,   2,   5,  NA,   7),
 mother = c( NA,  NA,   1,  NA,   3,   3,   6,  NA,   8),
 phenotype = c( NA, 105,  98, 101, 106,  93,  NA,  NA, 109),
 group = factor(c( NA,   1,   1,   2,   2,   2,  NA,  NA,   1))
)
tmp <- AnimalModelBUGSData(ped=example)
print(tmp)

}
\references{
Damgaard, L. H. Technical note: How to use Winbugs to
  draw inferences in animal models. \emph{Journal of Animal
  Science}, v. 85, n. 6, p. 1363-1368, 2007
}
\seealso{
\code{\link{writeJags}}
}
\author{
Thiago de Paula Oliveira, \email{thiago.oliveira@ed.ac.uk},
  Ivan Pocrnic, \email{ivan.pocrnic@roslin.ed.ac.uk}, Gregor Gorjanc,
  \email{gregor.gorjanc@gmail.com}
}
\keyword{BUGS}
\keyword{Mixed-Effects-Model}
\keyword{Pedigree}
