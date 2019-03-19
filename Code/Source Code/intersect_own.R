#' Fast Intersect Function R
#'
#' Faster version of base::intersect in R
#' @param a vector
#' @param b vector
#' @export


intersect_own <- function(a,b){
  c(a,b)[duplicated(c(a,b))]
}
#' Fast Intersect Function R
#'
#' Faster version of base::intersect in R
#' @param a vector
#' @param b vector
#' @export


intersect_own3 <- function(a,b){
  if(length(a)==0 || length(b)==0){
    return(NULL)
  }
  c(a,b)[duplicated(c(a,b))]
}

#' No-Fast Intersect Function R
#'
#' No-Faster version of base::intersect in R
#' @param a vector
#' @param b vector
#' @export

intersect_own2 <- function(a,b){
  activ_a <- 1
  activ_b <- 1
  activ_d <- 1
  la <- length(a)
  lb <- length(b)
  dups <- numeric(lb)
  while(activ_b<=lb && activ_a <=la){
    if(b[activ_b]<a[activ_a]){
      activ_b <- activ_b +1
    } else if(b[activ_b]>a[activ_a]){
      activ_a <- activ_a +1
    } else{
      dups[activ_d] <- b[activ_b]
      activ_b <- activ_b + 1
      activ_d <- activ_d + 1
      activ_a <- activ_a + 1
    }
  }
  if(activ_d==1){
    return(numeric(0))
  }
  return(dups[1:(activ_d-1)])
}
