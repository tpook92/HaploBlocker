#' Number of Indi
#'
#' Function to determine number of individuals in a blocklist
#' @param blocklist block-dataset
#' @examples
#' data(blocklist_ex_maze)
#' indi_calc(blocklist_ex_maze)
#' @export
#' @return Number of individuals considered in the haplotype library

indi_calc <- function(blocklist){
  indi <- 0
  for(index in 1:length(blocklist)){
    indi <- max(indi, blocklist[[index]][[6]])
  }

  return(indi)
}
