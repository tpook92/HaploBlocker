#' Determine in which block a single haplotype is in over the chromosom
#'
#' Determine in which block a single haplotype is in over the chromosom
#' @param data knot-dataset
#' @param indi number of haplotypes in the dataset
#' @param nwindow number of windows in the dataset
#' @export

current_block <- function(data, indi, nwindow){
  n <- length(data)
  blocknr <- numeric(nwindow)
  for(index in 1:n){
    check <- sum(data[[index]][[5]]==indi)
    if(check==1){
      blocknr[data[[index]][[1]]$window: data[[index]][[2]]$window] <- index
    }
  }
  return(blocknr)
}
