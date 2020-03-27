#' Determine in which node a single haplotype is in over the chromosom
#'
#' Determine in which node a single haplotype is in over the chromosom
#' @param data node-dataset
#' @param haplo_nr Which haplotype to check to node sequence
#' @param nwindow number of windows in the dataset
#' @export

current_block <- function(data, haplo_nr, nwindow){
  n <- length(data)
  blocknr <- numeric(nwindow)
  for(index in 1:n){
    check <- sum(data[[index]][[5]]==haplo_nr)
    if(check==1){
      blocknr[data[[index]][[1]]$window: data[[index]][[2]]$window] <- index
    }
  }
  return(blocknr)
}
