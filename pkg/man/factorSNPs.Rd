\name{factorSNPs}
\alias{factorSNPs}

\title{Factorization of the columns of a matrix}
\description{The functions determines which code in a given part of a SNP sequence is present in other individuals}

\usage{
factorSNPs(CM, start=1, end=attr(CM, "position")[2] + 1)
}
% hidden argument: check, verbose, mass.balancing, do.plot
\arguments{
  \item{CM}{a coded matrix of usually many columns and relatively few lines}
  \item{start,end}{start and end point in the SNP sequence}
}

\value{
  A vector with the factor value for each column (from 1 to the number
  of different columns).
  The attribute \code{"where.to.find"} gives the position where the
  different columns can be found for the first time.
  The attribute \code{"counts"} gives the number of columns found for
  each pattern.
}

\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}}

\seealso{
  \command{\link{codeSNPs}}, \command{\link{decodeSNPs}}
}

\examples{
%  library(CHaploBlocker)
require(RandomFieldsUtils)

fixcoding(0:1)
for (i in 1:10) {
  cat(i, "")
  animals <- sample(500, 1)
  snps <- sample(70, 1)
  v_snps <- sample(snps - 1, 1)
  start <- sample(snps - v_snps, 1)
  end <- sample(start:snps, 1)
  Print(animals, snps, start, end)
  
  M <- matrix(nc=animals, sample(c(0,1), animals * snps, replace = TRUE))
  M <- cbind(M, M)
  M <- M[, sample(ncol(M)), drop=FALSE]
  ##  print(M)
  CM <- codeSNPs(M)
  f <- factorSNPs(CM, start, end)
  ff <- unique(M[start:end, ,drop=FALSE], MARGIN = 2)
  ##
  stopifnot(ncol(ff) == length(attr(f, "where.to.find")))
  stopifnot(ncol(ff) == length(attr(f, "where.to.find")))
  stopifnot(all(ff == M[start:end, attr(f, "where.to.find")]))
}
}




\keyword{misc}
