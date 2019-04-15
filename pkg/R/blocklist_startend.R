#' Merge blocks
#'
#' Function to merge/remove/modify blocks
#' @param blocklist block-dataset
#' @param type length measure (default: "window" , alt: "snp", "bp")
#' @param first_block First block to consider in the computation (default: 1)
#' @export

blocklist_startend <- function(blocklist, type="snp", first_block=1){
  be <- matrix(0, nrow=(length(blocklist)-first_block+1), ncol=2)
  for(index in 1:(length(blocklist)-first_block+1)){
    be[index,1] <- blocklist[[index+first_block-1]][[2]][[type]]
    be[index,2] <- blocklist[[index+first_block-1]][[3]][[type]]
  }
  return(be)
}
