#' Help function for list length
#'
#' Function to derive the length of lists
#' @param list List
#' @examples
#' data(blocklist_ex_maze)
#' llength(blocklist_ex_maze)
#' @export
#' @return Length of each list element

llength <- function(list){
  len <- numeric(length(list))
  for(index in 1:length(list)){
    len[index] <- length(list[[index]])
  }
  return(len)
}
