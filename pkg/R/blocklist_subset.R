#' Merge blocks
#'
#' Function to merge/remove/modify blocks
#' @param blocklist block-dataset
#' @param subset Subset of haplotype to reduce the haplotype library to
#' @return haplotype library
#' @export

blocklist_subset <- function(blocklist, subset){

  for(index in length(blocklist):1){

    blocklist[[index]][[6]] <- intersect(blocklist[[index]][[6]], subset)
    blocklist[[index]][[5]] <- length(blocklist[[index]][[6]])

    if(blocklist[[index]][[5]]==0){
      blocklist[[index]] <- NULL
    } else{
      blocklist[[index]][[6]] <- which(duplicated(c(blocklist[[index]][[6]], subset))[-(1:blocklist[[index]][[5]])])
    }

  }

  return(blocklist)
}
