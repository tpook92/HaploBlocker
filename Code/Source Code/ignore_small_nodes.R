#' Ignore small nodes
#'
#' Internal Function remove windows which transition to less frequent haplotypes in adjecent windows
#' @param data node-dataset
#' @param indi number of haplotypes in the dataset
#' @param nwindow number of windows in the dataset
#' @param node_min minimum number of haplotypes per block (default: 5)
#' @param gap remove haplotypes in nodes adjacent to nodes with less than minimim_blocksize haplotypes in it (default: 10 windows)
#' @param intersect_func Used intersect-function (internally relevant for computation time)
#' @export

ignore_small_nodes <- function(data, indi , nwindow, node_min=5, gap=10, intersect_func=intersect){
  resttiere <- matrix(0, nrow=indi, ncol=nwindow)

  for(index in 1:length(data)){
    if(data[[index]][[3]]< node_min){
      resttiere[data[[index]][[5]], data[[index]][[1]]$window:data[[index]][[2]]$window] <- 1
      data[[index]][[3]] <- 0
      data[[index]][[5]] <- numeric(0)

      if(data[[index]][[2]]$window<nwindow){
        next1 <- unique(c(0,data[[index]][[6]][,1]))[-1]
        for(index2 in next1){
          data[[index2]] <- calculate_new_transition(data, index2, nwindow, intersect_func=intersect_func)
        }
      }


      if(data[[index]][[1]]$window>1){
        previous1 <- unique(c(0,data[[index]][[7]][,1]))[-1]
        for(index2 in previous1){
          data[[index2]] <- calculate_new_transition(data, index2, nwindow, intersect_func=intersect_func)
        }
      }

    }
  }




  # Additional 0er
  checker <- matrix(0, nrow=indi, ncol=nwindow)

  for(index in 1:indi){
    ones <- which(resttiere[index,]==1)
    if(sum(ones)>0){
      diff <- ones[-1] -ones[-length(ones)]
      reaction <- unique(c(0,((diff < gap) & (diff>1)) *(1:length(diff))))[-1]
      for(change in reaction){
        checker[index, (ones[change]+1):(ones[change+1]-1)] <- 1
      }
      if(ones[1]< gap && ones[1]>1){
        checker[index, 1:ones[1]] <- 1
      }
      if(ones[length(ones)]> (nwindow-gap+1) && ones[length(ones)]<nwindow){
        checker[index, (ones[length(ones)]+1):nwindow] <- 1
      }
    }
  }


  block1 <- matrix(0, nrow=indi, ncol=nwindow)
  for(index in 1:length(data)){
    block1[data[[index]][[5]], data[[index]][[1]]$window: data[[index]][[2]]$window] <- index
  }


  abc <- unique(c(0,block1 * checker))[-1]

  for(index in abc){
    removes <- which(checker[,data[[index]][[1]]$window]==1)
    data[[index]][[5]] <- intersect_func(data[[index]][[5]], (1:indi)[-removes])
    data[[index]][[3]] <- length(data[[index]][[5]])
    data[[index]] <- calculate_new_transition(data, index, nwindow, intersect_func=intersect_func)
  }

  data <-renaming_combi(data, nwindow)

  return(data)

}
