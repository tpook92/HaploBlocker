#' Calculate size of each block
#'
#' Function to calculate the size of each block
#' @param blocklist block-dataset
#' @param intersect_func Used intersect-function (internally relevant for computation time)
#' @param first_block First block to consider in the computation (default: 1)
#' @export

blocklist_size <- function(blocklist, intersect_func=intersect, first_block=1){
  be <- numeric(length(blocklist)-first_block+1)
  for(index in 1:(length(blocklist)-first_block+1)){
    be[index] <- blocklist[[index+first_block-1]][[5]]
  }
  return(be)
}

#' Calculate size of each block
#'
#' Function to calculate the size of each block
#' @param blocklist block-dataset
#' @param checker vector of haplotype nrs to consider
#' @param intersect_func Used intersect-function (internally relevant for computation time)
#' @export

blocklist_size_spes <- function(blocklist,checker,intersect_func=intersect){
  be <- numeric(length(blocklist))
  for(index in 1:length(blocklist)){
    be[index] <- length(intersect_func(blocklist[[index]][[6]], checker))
  }
  return(be)
}
