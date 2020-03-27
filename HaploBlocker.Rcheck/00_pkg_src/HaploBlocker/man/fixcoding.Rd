\name{fixcoding}
\alias{fixcoding}
\title{Determination of Coding}
\description{
  determines which kind of coding is given in the application
}
\usage{
fixcoding(values) 
}
% hidden argument: check, verbose, mass.balancing, do.plot
\arguments{
  \item{values}{character or integer.
    If a single integer is given, then coding for each SNP separately is
    performed. The modulus gives the number of maximal codes. If the
    value is positive and there are more different codes the most
    frequent ones are taken and the last value collects also all
    exceptions. If the value is negative no collection is performed
    and an error occurs if there are more codes.

    If \code{values} has more than one value then \code{values} gives
    the allowed codes (global coding). The last value can be
    \code{"ANY"}, then all the rest is collected here.

    Default is \code{fixcoding(values=-2)}.
  }
}
\details{  
  Per SNP only 256 different values are currently allowed.  
  If more SNPs values appear and \code{values} is a single positive
  value, only the first 256 SNPs are considered for coding.
  If \code{values} is negative and the 256 SNPs are exceeded, an error
  caused. From the maximal 256 SNPs the best (\code{values}-1) SNPs are
  coded if \code{values} is positive. 
}
\value{
  NULL
}
\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}}

\seealso{
  \command{\link{codeSNPs}}, \command{\link{decodeSNPs}}
}

\examples{
fixcoding(0:1)

animals <- 8
snps <- 25
(M <- matrix(nc=animals, sample(0:1, animals*snps, replace=TRUE)))
(CM <- codeSNPs(M))
(decCM <- decodeSNPs(CM))
confirm(M, decCM)

}

\keyword{misc}
