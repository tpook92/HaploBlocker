#' Comparison of blocklists
#'
#' Function compair different blocklists
#' @param blocklist1 blocklist created by blocklist_calculcation algorithm
#' @param blocklist2 blocklist created by blocklist_calculcation algorithm
#' @param type length measure (default: "window" , alt: "snp", "bp")
#' @param shift switch the positions of the second blocklist (default: 0)
#' @param turn.around If TRUE turn the results of blocklist2 around
#' @param indi number of haplotypes in the dataset
#' @param compair_region Region to compair (default: Every overlapping SNP)
#' @param intersect_func Used intersect-function (internally relevant for computation time)
#' @export

blocklist_comparison<- function(blocklist1, blocklist2, indi=NULL, type="snp", compair_region=NULL, shift=0, turn.around=FALSE, intersect_func=intersect){

  if(length(indi)==0){
    indi <- indi_calc(blocklist1)
  }
  se1 <- blocklist_startend(blocklist1, type=type)
  se2 <- blocklist_startend(blocklist2, type=type)
  if(turn.around==TRUE){

    se2 <- max(se2[,2]) - se2 + 1
    se2[,1:2] <- se2[,2:1]
  }
  se2 <- se2 + shift
  if(length(compair_region)==0){
    compair_region <- max(min(se1[,1]), min(se2[,1])):min(max(se1[,2]), max(se2[,2]))
  }
  sim_score <- matrix(0, ncol=length(compair_region), nrow=indi)
  min_s <- min(compair_region)
  sim <- matrix(0, nrow=length(blocklist1), ncol=length(blocklist2))
  for(index1 in 1:length(blocklist1)){
    for(index2 in 1:length(blocklist2)){
      inters <- intersect_func(blocklist1[[index1]][[6]], blocklist2[[index2]][[6]])
      sim[index1,index2] <- length(inters) / length(unique(c(blocklist1[[index1]][[6]], blocklist2[[index2]][[6]])))
      if(max(se1[index1,1], se2[index2,1])<=min(se1[index1,2], se2[index2,2])){
        overlap <- max(se1[index1,1], se2[index2,1]):min(se1[index1,2], se2[index2,2])
      } else{
        overlap <- NULL
      }
      overlap <- intersect_func(overlap, compair_region)
      old <- sim_score[inters,overlap-min_s+1]
      old[old<sim[index1,index2]] <- sim[index1,index2]
      sim_score[inters,overlap-min_s+1] <- old
    }
  }
  t1 <- coverage_test(blocklist1, indi, type=type)[compair_region,]
  t2 <- coverage_test(blocklist2, indi, type=type)
  if(turn.around==TRUE){
    t2 <- t2[nrow(t2):1,]
  } else{
    t1 <- coverage_test(blocklist1, indi, type=type)[compair_region+shift,]
  }
  t2 <- t2[compair_region,]
  t_sum <- t(t1+t2)
  sim_score[t_sum==0] <- 1
  return(sim_score)
}
