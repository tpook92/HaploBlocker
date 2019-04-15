#' Number of Indi
#'
#' Function to determine number of individuals in a blocklist
#' @param blocklist block-dataset
#' @export

indi_calc <- function(blocklist){
  indi <- 0
  for(index in 1:length(blocklist)){
    indi <- max(indi, blocklist[[index]][[6]])
  }

  return(indi)
}
