#' Function to plot a blocklist II
#'
#' Function to plot a blocklist. Width of blocks is adjusted according to nr. of haplotypes in the block
#' @param blocklist blocklist created by blocklist_calculcation algorithm
#' @param cutoff2 minimum number of blocks to start/end to mark a recombination hotspot (default:5)
#' @param bound_weighted weighted blocks in the detection of recombination hotspots according to size (default: TRUE)
#' @param type length measure (default: "window" , alt: "snp", "bp")
#' @export

blocklist_plot_xsize <- function(blocklist, cutoff2 = 5, bound_weighted=TRUE, type="snp" ){
  size <- blocklist_size(blocklist)
  se <- blocklist_startend(blocklist, type=type)
  nwindow <- max(blocklist_startend(blocklist, type="window"))

  bpstart <- se[,1]
  bpend <- se[,2]

  switch <- sort(bpstart, decreasing=FALSE, index.return=TRUE)$ix
  start <- bpstart[switch]
  end <- bpend[switch]
  size <- size[switch]
  plot(0,-100, xlim=c(1,sum(size)), ylim=c(0, max(end)), ylab=type, xlab="Indi per Block")
  current_value <- 1
  for(index in 1:length(start)){
    polygon(c(current_value,current_value+size[index], current_value+ size[index],current_value), c(start[index], start[index], end[index], end[index]))
    current_value <- current_value + size[index]
  }

  if(bound_weighted==TRUE){
    tempest <- numeric(sum(size)*2)
    sofar <- 0
    for(index in 1:length(start)){
      tempest[1:size[index]+sofar] <- start[index]
      sofar <- sofar + size[index]
    }
    for(index in 1:length(start)){
      tempest[1:size[index]+sofar] <- end[index]
      sofar <- sofar + size[index]
    }
    xy <- density(tempest,bw=5, from=1, to=max(se), n=nwindow)

    candidate <- unique(c(0,(c(0,xy$y[-nwindow])<c(xy$y)) * (c(xy$y)>c(xy$y[-1],0)) * (xy$y > 1/length(start)/2/5/2*cutoff2) * 1:nwindow))[-1]

    abline(h=candidate)
  } else{
    xy <- density(c(start,end),bw=5, from=1, to=nwindow, n=nwindow)

    candidate <- unique(c(0,(c(0,xy$y[-nwindow])<c(xy$y)) * (c(xy$y)>c(xy$y[-1],0)) * (xy$y > 1/length(start)/2/5/2*cutoff2) * 1:nwindow))[-1]

    abline(h=candidate)
  }
}
