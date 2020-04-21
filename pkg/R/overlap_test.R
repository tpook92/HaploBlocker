#' Overlap Test
#'
#' Function to calculate the overlap of the blocks with another dataset
#' @param blocklist block-dataset
#' @param dhm dataset to check for overlap with the blocklist
#' @param max maximum number in each cell (overlapping blocks; default: 1)
#' @param min_similarity minimum rate of the same SNPs to be added to the block (default: 0.99)
#' data(ex_maze)
#' data(blocklist_ex_maze)
#' overlap <- overlap_test(blocklist_ex_maze, ex_maze)
#' @export
#' @return Overlap matrix

overlap_test <- function(blocklist, dhm , max=1, min_similarity=0.99){
  overlap <- matrix(0, nrow=nrow(dhm), ncol=ncol(dhm))
  for(index in 1:length(blocklist)){
    seq <- blocklist[[index]][[7]]$snp
    start <- blocklist[[index]][[2]]$snp
    end <- blocklist[[index]][[3]]$snp

    sim <-colSums(dhm[start:end,]== seq)
    included <- which(sim >= (min_similarity*(end-start+1)))
    if(length(included)>0){
      overlap[start:end, included] <- overlap[start:end, included] +1
    }
  }
  overlap[overlap>max] <- max

  return(overlap)
}
