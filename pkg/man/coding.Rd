\name{coding}
\alias{coding}
\alias{codeSNPs}
\alias{decodeSNPs}
\title{Coding}
\description{
  These function code and decode genetic matrices according to the code
  fixed by  \command{\link{fixcoding}}
}
\usage{
codeSNPs(M, redo=is.character(M), SNPxINDIVID=TRUE)

decodeSNPs(CM)
}
% hidden argument: check, verbose, mass.balancing, do.plot
\arguments{
  \item{M}{Genetik matrix or vector to be coded}
  \item{redo}{logical. If \code{TRUE} the current matrix is used to
    define further (exceptional) codes, if new code appear in \code{M}
    that have not been defined or observed yet.}
  \item{SNPxINDIVID}{currently disabled! It will allow to read in also
    transposed matrices}
  \item{CM}{coded matrix or vector}
}
\value{
  \command{codeSNPs} returns a coded matrix or a vector.
  \command{decodeSNPs} returns a matrix or a vector in human readable form. 
}

\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}}

\seealso{
  \command{\link{fixcoding}}
}

\examples{
animals <- 8
snps <- 25
ACGT <- c("A", "C", "G", "T")
fixcoding(c(ACGT, "ANY"))
(M <- matrix(nc=animals, sample(c(ACGT, "-", "+"), animals*snps, replace=TRUE)))
(CM <- codeSNPs(M))
(decCM <- decodeSNPs(CM))
stopifnot(all(M == decCM || ((M == "-" || M == "+") && decCM == "@") ))

%(CM <- codeSNPs(M, SNP=FALSE))
%(decCM <- decodeSNPs(CM))
%stopifnot(all(M == decCM || ((M == "-" || M == "+") && decCM == "@") ))
}


\keyword{models}
