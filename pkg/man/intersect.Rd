\name{intersect_c}
\alias{intersect_c}
\title{Intersection of sets}
\description{
  A fast implementation of R command \command{\link[base]{intersect}}
  for interger valued, sorted vectors
}
\usage{
intersect_c(a, b)
}
% hidden argument: check, verbose, mass.balancing, do.plot
\arguments{
  \item{a,b}{integer valued, sorted vectors representing sets}
}

\value{
  sorted vector of integer values representing the intersection of the
  sets \code{a} and \code{b}.
}

\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}}

\seealso{
  \command{\link[base]{intersect}}
}

\examples{
(a <- as.integer(sort(sample(20, 6))))
(b <- as.integer(sort(sample(20,7))))
confirm(intersect_c(a,b), base::intersect(a,b))

a <- as.integer(sort(sample(2000, 600)))
b <- as.integer(sort(sample(2000,700)))
confirm(intersect_c(a,b), base::intersect(a,b))
print(system.time(for (i in 1:10^5) base::intersect(a,b)))
print(system.time(for (i in 1:10^5) intersect_c(a,b))) ## factor 10 faster

}

\keyword{misc}
