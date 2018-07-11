#' nodes size
#'
#' Function to calculate the number of haplotypes in each node
#' @param data node-dataset
#' @export

nodes_size <- function(data){
  n <- length(data)
  size <- numeric(n)
  for(index in 1:n){
    size[index] <- data[[index]][[3]]
  }
  return(size)
}
