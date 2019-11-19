#' Help function for list length
#'
#' Function to derive the length of lists
#' @param list List
#' @export

llength <- function(list){
  len <- numeric(length(list))
  for(index in 1:length(list)){
    len[index] <- length(list[[index]])
  }
  return(len)
}
