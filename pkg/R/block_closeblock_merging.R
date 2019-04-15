#' Merging neighbouring blocks
#'
#' Function to merging neighbouring block with similar haplotypes/allels
#' @param blocklist block-dataset
#' @param blockinfo List with all relevant information to each window seperatly
#' @param indi number of haplotypes in the dataset
#' @param nwindow number of windows in the dataset
#' @param intersect_func Used intersect-function (internally relevant for computation time)
#' @param max_diff_l maximum number of windows with different haplotypes inbetween (default: 1)
#' @param max_diff_i maximum number of individuals in only one of the two blocks (default: 1)
#' @param dataset dataset which variant nr. for each window
#' @export

block_closeblock_merging <- function(blocklist, blockinfo,  indi, nwindow, max_diff_l, max_diff_i, intersect_func=intersect, dataset){
  same <- matrix(0, ncol=length(blocklist), nrow=length(blocklist))

  for(index in 1:length(blocklist)){
    for(index2 in index:length(blocklist)){
      same[index,index2] <- same[index2,index] <- length(intersect_func(blocklist[[index]][[6]], blocklist[[index2]][[6]]))
    }
  }
  merges <- (same>=diag(same)-max_diff_i) * t(same>=diag(same)-max_diff_i)
  a <- which(merges==1, arr.ind=TRUE)
  a <- a[a[,1]!=a[,2],]
  a <- a[a[,1]>a[,2],,drop=FALSE]
  a <- rbind(rbind(a,NULL)[,2:1],NULL)

  if(nrow(a)>0){
    for(index in 1:nrow(a)){
      vor <- a[index,1]
      nach <- a[index,2]
      cluster <- blocklist[[index]][[12]]
      if(vor>0){
        vor_ende <- blocklist[[vor]][[3]]$window
        nach_start <- blocklist[[nach]][[2]]$window

        mittelwerte <- dataset[[cluster]][intersect_func(blocklist[[vor]][[6]], blocklist[[nach]][[6]]), -c(1:vor_ende, nach_start:nwindow), drop=FALSE]
        variants <- numeric(ncol(mittelwerte))
        common <- numeric(ncol(mittelwerte))
        if(length(variants)>0){
          for(index2 in 1:length(variants)){
            variants[index2] <- length(unique(mittelwerte[,index2]))
            if(variants[index2]==1){
              common[index2] <- mittelwerte[1,index2]
            } else{
              versions <- unique(mittelwerte[,index2])
              counts <- numeric(length(versions))
              for(index3 in 1:length(versions)){
                counts[index3] <- sum(mittelwerte[,index2]==versions[index3])
              }
              common[index2] <- versions[which.max(counts)[1]]
            }
          }
        }

        if(sum(variants>1)<= max_diff_l){
          blocklist[[vor]][[1]] <- sort(unique(blocklist[[vor]][[1]], blocklist[[nach]][[1]]))
          prev <- blocklist[[vor]][[3]]$window - blocklist[[vor]][[2]]$window
          blocklist[[vor]][[2]]$window <- min(blocklist[[vor]][[2]]$window, blocklist[[nach]][[2]]$window)
          blocklist[[vor]][[3]]$window <- max(blocklist[[vor]][[3]]$window, blocklist[[nach]][[3]]$window)
          blocklist[[vor]][[2]]$snp <- min(blocklist[[vor]][[2]]$snp, blocklist[[nach]][[2]]$snp)
          blocklist[[vor]][[3]]$snp <- max(blocklist[[vor]][[3]]$snp, blocklist[[nach]][[3]]$snp)
          blocklist[[vor]][[2]]$bp <- min(blocklist[[vor]][[2]]$bp, blocklist[[nach]][[2]]$bp)
          blocklist[[vor]][[3]]$bp <- max(blocklist[[vor]][[3]]$bp, blocklist[[nach]][[3]]$bp)
          blocklist[[vor]][[6]] <- intersect_func(blocklist[[vor]][[6]], blocklist[[nach]][[6]])
          blocklist[[vor]][[5]] <- length(blocklist[[vor]][[6]])
          if(prev+length(common) == (blocklist[[vor]][[3]]$window-blocklist[[vor]][[2]]$window)){
            blocklist[[vor]][[4]] <- c(blocklist[[vor]][[4]], common)

          } else{
            blocklist[[vor]][[4]] <- c(blocklist[[vor]][[4]], common, blocklist[[nach]][[4]][(length(blocklist[[nach]][[4]])- (blocklist[[vor]][[3]]$window-blocklist[[vor]][[2]]$window-length(blocklist[[vor]][[4]])-length(common))):length(blocklist[[nach]][[4]])])
          }
          change <- which(a[,1]==a[index,2])
          new_diff <- which(variants>1) + vor_ende

          if(length(blocklist[[vor]])>=8 && length(blocklist[[nach]])>=8){
            blocklist[[vor]][[8]] <- sort(c(blocklist[[vor]][[8]], blocklist[[nach]][[8]], new_diff))
          } else if(length(blocklist[[vor]])>=8){
            blocklist[[vor]][[8]] <- sort(c(blocklist[[vor]][[8]], new_diff))
          } else if(length(blocklist[[nach]])>=8){
            blocklist[[vor]][[8]] <- sort(c(blocklist[[nach]][[8]], new_diff))
          } else{
            blocklist[[vor]][[8]] <- sort(new_diff)
          }
          a[change,1] <- 0
        } else{
          a[index,2] <- 0
        }
      }
    }
  }
  removes <- unique(c(0,a[,2]))[-1]
  for(rindex in sort(removes, decreasing = TRUE)){
    blocklist[[rindex]] <- NULL
  }
  return(blocklist)
}

