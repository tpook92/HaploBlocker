#' Old Major Version Identifier
#'
#' Old internal function to identify the major version of each haplotype in each window
#' @param blockinfo List with all relevant information to each window seperatly
#' @param dhm haploid SNP-dataset
#' @param window_size size of each window in the algorithm (default: 20)
#' @export

blockinfo_major <- function(blockinfo, dhm, window_size){
  for(index in 1:length(blockinfo)){
    blockinfo[[index]][[6]] <- list()
    for(index2 in 1:length(blockinfo[[index]][[4]])){
      count <- numeric(nrow(blockinfo[[index]][[4]][[index2]]))
      for(index3 in 1:length(count)){
        count[index3] <- sum(dhm[(window_size*(index-1)+1):(window_size*index), blockinfo[[index]][[5]][[index2]]] == blockinfo[[index]][[4]][[index2]][index3,])
      }
      blockinfo[[index]][[6]][[index2]] <- blockinfo[[index]][[4]][[index2]][which.max(count),]

    }
  }
  return(blockinfo)
}
