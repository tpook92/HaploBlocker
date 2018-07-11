#' Utility Unique Matrix
#'
#' Internal Function for easy unique function on matrixes
#' @param mat Matrix
#' @export


uniquem <- function(mat){
  vars <- unique(mat)
  counts <- numeric(nrow(vars))
  for(index in 1:nrow(vars)){
    counts[index] <- sum(colSums(vars[index,] == t(mat))==ncol(vars))
  }
  return(list(vars, counts))
}
