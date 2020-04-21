#' Utility Unique Matrix
#'
#' Internal Function for easy unique function on matrixes
#' @param mat Matrix
#' @examples
#' dataset <- rbind(c(1,1,1), c(1,2,1), c(1,2,1), c(2,2,2))
#' uniquem(dataset)
#' @export
#' @return [[1]] present variants [[2]] frequency


uniquem <- function(mat){
  vars <- unique(mat)
  counts <- numeric(nrow(vars))
  for(index in 1:nrow(vars)){
    counts[index] <- sum(colSums(vars[index,] == t(mat))==ncol(vars))
  }
  return(list(vars, counts))
}
