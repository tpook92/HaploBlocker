#' knots size
#'
#' Function to calculate the number of haplotypes in each knot
#' @param data knot-dataset
#' @export

knoten_size <- function(data){
  n <- length(data)
  size <- numeric(n)
  for(index in 1:n){
    size[index] <- data[[index]][[3]]
  }
  return(size)
}
