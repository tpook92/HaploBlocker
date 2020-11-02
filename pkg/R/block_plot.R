#' Plotfunction for the block library
#'
#' Function to plot the location and size of block and the resulting coverage
#' @param blocklist block-dataset
#' @param type length measure (default: "window" , alt: "snp")
#' @param indi number of haplotypes in the dataset
#' @param bw bandwidth for the smoothing of the coverage (default: bw=1 - no smoothing)
#' @examples
#' data(blocklist_ex_maze)
#' block_plot(blocklist_ex_maze)
#' @export
#' @return Visualization of local coverages/block variation

block_plot <- function(blocklist,indi=NULL, type="snp", bw=1){

  if(length(indi)==0){
    indi <- indi_calc(blocklist)
  }

  se <- blocklist_startend(blocklist, type=type)
  size <- blocklist_size(blocklist)/indi
  xlab <- "SNP"
  if(type=="bp"){
    xlab = "base pair"
  }
  if(type=="window"){
    xlab = "window"
  }


  plot(0,-100, xlim=c(0, max(se[,2])), ylim=c(0, 1), ylab="Size/Coverage", xlab=xlab, main="Location of the blocks")
  for(index in 1:length(blocklist)){
    lines(se[index,c(1,2,2,1,1)],c(0,0,size[index], size[index],0))
  }
  if(type!="bp"){
    cov <- rowMeans(coverage_test(blocklist, indi = indi, type=type))
  } else{
    cov <- rowMeans(coverage_test(blocklist, indi = indi, type="snp"))
  }

  lines(ksmooth(1:length(cov), cov, bandwidth = bw), col="red")

}
