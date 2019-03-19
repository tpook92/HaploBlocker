#' Plotfunction for the block library
#'
#' Function to plot the location and size of block and the resulting coverage
#' @param blocklist block-dataset
#' @param type length measure (default: "window" , alt: "snp", "bp")
#' @param indi number of haplotypes in the dataset
#' @param bw bandwidth for the smoothing of the coverage (default: bw=1 - no smoothing)
#' @export

block_plot <- function(blocklist,indi=313, type="snp", bw=1){
  se <- blocklist_startend(blocklist, type=type)
  size <- blocklist_size(blocklist)/indi
  plot(0,-100, xlim=c(0, max(se[,2])), ylim=c(0, 1), ylab="Size/Coverage", xlab=type, main="Location of the blocks")
  for(index in 1:length(blocklist)){
    lines(se[index,c(1,2,2,1,1)],c(0,0,size[index], size[index],0))
  }
  cov <- rowMeans(coverage_test(blocklist, indi = indi, type=type))
  lines(ksmooth(1:length(cov), cov, bandwidth = bw), col="red")

}
