#' Generate Blockmatrix
#'
#' Function to create a block-matrix
#' @param blocklist blocklist created by blocklist_calculcation algorithm
#' @param indi number of haplotypes in the dataset
#' @export

generate_blockmatrix <- function(blocklist, indi){
  blockmatrix <- matrix(0, ncol=indi, nrow=length(blocklist))
  for(index in 1:length(blocklist)){
    blockmatrix[index,blocklist[[index]][[6]]] <- 1
  }
  return(blockmatrix)
}
