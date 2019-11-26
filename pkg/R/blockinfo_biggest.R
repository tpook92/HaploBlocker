#' Compute biggest block for each region
#'
#' Internal function to count positions where each block is major block
#' @param blocklist block-dataset
#' @param indi number of haplotypes in the dataset
#' @param nwindow number of windows in the dataset
#' @param type length measure (default: "window" , alt: "snp", "bp")
#' @param min_majorblock minimum of positions in the dataset a block is the biggest covering (default: 50)
#' @param weighting_length Weighting factor for length to determine major block (default: 1)
#' @param weighting_size Weighting factor for number of haplotypes in block to determine major block (default: 1)
#' @param recalculate_biggest Set to FALSE to only calculate the number of major positions for those blocks that could be removed in each iteration (number of major blocks is only increasing when removing other blocks)
#' @param window_size size of each window in the algorithm (default: 20)
#' @export

blockinfo_biggest <- function(blocklist, nwindow=NULL, indi=NULL, type="window", min_majorblock=5000, weighting_length=1, weighting_size=1,
                              recalculate_biggest=TRUE, window_size, deletion_count=FALSE, present_data=NULL){
  if(length(unique(window_size))!=1){
    type <- "snp"
  } else{
    min_majorblock <- min_majorblock / window_size[1]
  }

  if(length(indi)==0){
    indi <- indi_calc(blocklist)
  }
  if(length(nwindow)==0){
    helper <- max(blocklist_startend(blocklist, type="snp"))
    nwindow <- ceiling(helper/window_size)
  }

  se <- blocklist_startend(blocklist, type=type)
  sel <- se[,2]-se[,1]+1

  if(deletion_count){
    for(index in 1:length(blocklist)){
      sel[index] <- sum(blocklist[[index]][[7]]$snp!=0) / window_size[1]
    }
  }

  size <- blocklist_size(blocklist)
  major_rating <- sel ^ weighting_length * size ^ weighting_size
  order <- sort(major_rating, index.return=TRUE)$ix

  needed_size <- max(se[,2])
  bdataset <- matrix(0, nrow=indi, ncol=needed_size)

  for(index in order){
    if(type=="window"){
      bdataset[blocklist[[index]][[6]], blocklist[[index]][[2]]$window:blocklist[[index]][[3]]$window] <- index
    } else{
      bdataset[blocklist[[index]][[6]], blocklist[[index]][[2]]$snp:blocklist[[index]][[3]]$snp] <- index
    }

  }

  count <- numeric(length(blocklist))
  for(index in 1:length(count)){
    if(recalculate_biggest || length(blocklist[[index]])<11 || length(blocklist[[index]][[11]])==0 || blocklist[[index]][[11]] < min_majorblock){
      if(type=="window"){
        if(deletion_count){
          count[index] <- sum((bdataset[blocklist[[index]][[6]],blocklist[[index]][[2]]$window:blocklist[[index]][[3]]$window]==index)*
                                present_data[blocklist[[index]][[6]], blocklist[[index]][[2]]$window:blocklist[[index]][[3]]$window])
        } else{
          count[index] <- sum(bdataset[blocklist[[index]][[6]],blocklist[[index]][[2]]$window:blocklist[[index]][[3]]$window]==index)
        }

      } else{
        count[index] <- sum(bdataset[blocklist[[index]][[6]],blocklist[[index]][[2]]$snp:blocklist[[index]][[3]]$snp]==index)
      }

      blocklist[[index]][[11]] <- count[index]
    } else{
      count[index] <- blocklist[[index]][[11]]
    }

  }
  remove_count <- sort(which(count<min_majorblock), decreasing=TRUE)

  for(index in remove_count){
    blocklist[[index]] <- NULL
  }
  return(blocklist)
}
