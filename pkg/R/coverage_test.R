#' Coverage Test
#'
#' Function to calculate the coverage of the blocks over the whole dataset
#' @param blocklist block-dataset
#' @param indi number of haplotypes in the dataset
#' @param type length measure (default: "window" , alt: "snp", "bp")
#' @param max maximum number in each cell (overlapping blocks; default: 1)
#' @examples
#' data(blocklist_ex_maze)
#' t <- coverage_test(blocklist_ex_maze)
#' @export
#' @return coverage matrix

coverage_test <- function(blocklist, indi=NULL, type="snp", max=1){
  if(length(indi)==0){
    indi <- indi_calc(blocklist)
  }

  se <- blocklist_startend(blocklist, type=type )
  nwindow <- max(se[,2])
  coverage.test <- matrix(0, ncol=indi, nrow=nwindow)
  for(index in 1:length(blocklist)){
    coverage.test[blocklist[[index]][[2]][[type]]:blocklist[[index]][[3]][[type]], blocklist[[index]][[6]]] <- coverage.test[blocklist[[index]][[2]][[type]]:blocklist[[index]][[3]][[type]], blocklist[[index]][[6]]] +1
  }
  coverage.test[coverage.test>max] <- max

  return(coverage.test)
}

