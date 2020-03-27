\name{colSumsEqualSNPs}
\alias{colSumsEqualSNPs}
\title{Columwise comparison of two genetic matrices}
\description{
  For each column of the matrix \code{CM} and for each column of the
  matrix {V}. The number of SNPs in \code{CM} that are equal to those of
  \code{V} are calculated.
}
\usage{
colSumsEqualSNPs(CM, start=1, V, select = NULL)
}
% hidden argument: check, verbose, mass.balancing, do.plot
\arguments{
  \item{CM}{a (coded) genetik matrix or vector}
  \item{start}{socket row for \code{CM} and \code{V} where comparison
    starts.
    NOTE: If \code{V} (or \code{CM}) is given as a numerical
    vector, then it is expected that the \code{V[1]} (or \code{CM[1]}) is the value of the
    start position.
    If \code{V} is a coded vector then \code{V} must have
    been coding from position 0 on.
  }
  \item{V}{matrix or vector.
    \code{V} may have \code{NA}s at the beginning or at the end.
    If \code{V} is not given, \code{CM} is used instead}
  \item{select}{\code{select} gives the selected columns within
    \code{CM}. if not given all columns are of \code{CM} are selected.}
}
\details{
  if \code{CM} is identical to \code{V} then a symmetric matrix is
  returned by assumption (acceleration by factor 2 of the calculations).
}
\value{
  If \code{V} is a vector then a vector of length \code{nrow(CM)} is
  returned.
  Otherwise a matrix of size \code{nrow(CM)}\eqn{\times}\code{nrow(CM)}
  is returned.
}
\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}}

\seealso{
  \command{\link{codeSNPs}}, \command{\link{decodeSNPs}}
}

\examples{
fixcoding(0:1)
require(RandomFieldsUtils)
for (i in 1:10) {
  animals <- sample(100, 1)
  snps <- sample(30, 1)
  v_snps <- sample(snps, 1)
  start <- sample(snps - v_snps + 1, 1)

  M <- matrix(nc=animals, sample(c(0,1), animals * snps, replace = TRUE))
  V <- sample(c(0,1), v_snps, replace = TRUE)
  Vext <- c(rep(NA, start-1), V, rep(NA, snps - v_snps - start + 1))
  stopifnot(length(Vext) == nrow(M))
  CM <- codeSNPs(M)
  e <- colSumsEqualSNPs(CM, start, V)
  Print(animals, snps, v_snps, start, e)
  confirm(e, colSums(M == Vext, na.rm=TRUE))
}

}


\keyword{misc}
