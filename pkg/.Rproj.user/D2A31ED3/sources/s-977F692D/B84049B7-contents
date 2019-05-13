#' Which block
#'
#' Function to derivde in which block a giving haplotype is
#' @param blocklist block-dataset
#' @param nr number of haplotype to analyse
#' @export

which.block <- function(blocklist, nr){
  include <- numeric(length(blocklist))
  for(index in 1:length(blocklist)){
    if(sum(blocklist[[index]][[6]]==nr)>0){
      include[index] <-1
    }

  }
  return(include)
}

#' Which indi
#'
#' Function to derivde in which haplotypes are locally in the same block as a giving haplotype
#' @param blocklist block-dataset
#' @param nr number of haplotype to analyse
#' @param start Lower boundary of the considered intervall
#' @param end Upper boundary of the considered intervall
#' @export

which.indi <- function(blocklist, nr, start, end){
  indis <- NULL
  for(index in 1:length(blocklist)){
    if(blocklist[[index]][[2]]$bp<=start && blocklist[[index]][[3]]$bp>=end && sum(blocklist[[index]][[6]]==nr)>0){
      indis <- sort(unique(c(indis, blocklist[[index]][[6]])))
    }
  }
  return(indis)
}
