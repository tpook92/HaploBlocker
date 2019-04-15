#' Extract size of nodes
#'
#' Function to extract the size of the nodes of the node-dataset
#' @param data node-dataset
#' @export

data_size <- function(data){
  size <- numeric(length(data))
  for(index in 1:length(data)){
    size[index] <- data[[index]]$n_haplo
  }
  return(size)
}
