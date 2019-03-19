#' Coverage Test
#'
#' Function to calculate the coverage of the blocks over the whole dataset
#' @param blocklist block-dataset
#' @param indi number of haplotypes in the dataset
#' @param type length measure (default: "window" , alt: "snp", "bp")
#' @param max maximum number in each cell (overlapping blocks; default: 1)
#' @export

coverage_test <- function(blocklist, indi, type="snp", max=1){
  se <- blocklist_startend(blocklist, type=type )
  nwindow <- max(se[,2])
  coverage.test <- matrix(0, ncol=indi, nrow=nwindow)
  for(index in 1:length(blocklist)){
    coverage.test[blocklist[[index]][[2]][[type]]:blocklist[[index]][[3]][[type]], blocklist[[index]][[6]]] <- coverage.test[blocklist[[index]][[2]][[type]]:blocklist[[index]][[3]][[type]], blocklist[[index]][[6]]] +1
  }
  coverage.test[coverage.test>max] <- max

  return(coverage.test)
}


#blocktest <- numeric(length(blocklist))

#for(index in 1:length(blocklist)){
#  if(sum(blocklist[[index]][[6]]==9)&& blocklist[[index]][[2]]$window<=760&& blocklist[[index]][[3]]$window>760){
##    blocktest[index] <- 1
#  }
#}
