#' Sort Blocklist
#'
#' Function to sort the blocks based on starting-position
#' @param blocklist block-dataset
#' @param node_min minimum number of haplotypes per block (default: 5)
#' @return haplotype library


blocklist_reorder <- function(blocklist, node_min=5){
  size <- blocklist_size(blocklist)
  removes <- which(size<node_min)
  if(length(removes)>0){
    for(index in sort(removes, decreasing=TRUE)){
      blocklist[[index]] <- NULL
    }
  }
  se <- blocklist_startend(blocklist, type="snp")
  new.order <- sort(se[,1],index.return=TRUE)$ix
  new_blocklist <- list()
  for(index in 1:length(blocklist)){
    new_blocklist[[index]] <- blocklist[[new.order[index]]]
  }
  return(new_blocklist)
}
