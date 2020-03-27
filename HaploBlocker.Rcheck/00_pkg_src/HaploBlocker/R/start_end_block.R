#' node-boundaries
#'
#' Function to calculate start and end points of a node dataset
#' @param data node-dataset
#' @export


start_end_block <- function(data){
  n <- length(data)
  start <- numeric(n)
  end <- numeric(n)
  for(index in 1:n){
    start[index] <- data[[index]][[1]]$window
    end[index] <- data[[index]][[2]]$window
  }
  return(cbind(start,end))
}
