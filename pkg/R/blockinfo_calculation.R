#' Mainfunction to calculate haplotype blocks
#'
#' Function to generate haplotype blocks from haploid data
#' @param dhm haploid SNP-dataset
#' @param c_dhm Bit-wise coded SNP-dataset
#' @param window_size size of each window in the algorithm (default: 20)
#' @param merging_error number of allowed errors per block (default: 1)
#' @param window_sequence sequence of predefined windows (default: NULL ;per row: start$snp, end$snp, length, length - merging_error, start$bp, end$bp)
#' @param bp_map vector of positions for each SNP in bp (default: NULL - all 0)
#' @param window_anchor_gens matrix to constructed window_sequence base on start/end points in bp (e.g. gen regions, per row: start, end)
#' @param at_least_one If TRUE no allowed merging errors in windows of size 1
#' @param blockinfo_mode Structure of the groups in step I (default: 0 - Common haplos as major variants, 1- minimum number of groups)
#' @param max_groups Maximum number of groups per window (adaptive window-size, default: 0 - use fixed window size)
#' @param verbose Set to FALSE to not display any prints
#' @return [[1]] window-based haplotype information [[2]] window sequence

blockinfo_calculation <- function(dhm, window_sequence=NULL, window_size=20, merging_error=1, bp_map=NULL,
                                  window_anchor_gens=NULL, at_least_one=TRUE, max_groups=0, blockinfo_mode=0,
                                  c_dhm=NULL, verbose=TRUE){

  if(length(window_anchor_gens)>0 && length(window_sequence)==0 && length(bp_map)>0){
    window_sequence <- matrix(0, ncol=2, nrow=nrow(window_anchor_gens)*2+1)
    start <- 1
    for(index in 1:nrow(window_anchor_gens)){
      end <- sum(bp_map<window_anchor_gens[index,1])
      window_sequence[index*2-1,] <- c(start, end)
      end2 <- sum(bp_map<=window_anchor_gens[index,2])
      window_sequence[index*2,] <- c(end+1, end2)
      start <- end2+1
    }
    window_sequence[index*2+1,] <- c(start, length(bp_map))
    removes <-(window_sequence[,1]>window_sequence[,2])*(1:nrow(window_sequence))
    window_sequence <- window_sequence[-removes,]
  }

  if(max_groups==0 && length(window_sequence)==0){
    window_sequence <- cbind(0:(nrow(dhm)/window_size-1)*window_size+1, 1:(nrow(dhm)/window_size)*window_size)

    if( nrow(dhm)/window_size> nrow(window_sequence)){
      window_sequence <- rbind(window_sequence, c(window_size * nrow(window_sequence) + 1, nrow(dhm)))
    }
  }

  if(max_groups==0 && ncol(window_sequence)==2){
    window_sequence <- cbind(window_sequence, window_sequence[,2]-window_sequence[,1]+1, window_sequence[,2]-window_sequence[,1]-merging_error+1)
  }
  if(max_groups==0 && ncol(window_sequence)==3){
    window_sequence <- cbind(window_sequence[,1:2], window_sequence[,2]-window_sequence[,1]+1, window_sequence[,2]-window_sequence[,1]+1-window_sequence[,3])
  }
  if(max_groups==0 && at_least_one== TRUE && sum(window_sequence[,4]==0)>0){
    window_sequence[window_sequence[,4]==0, 4] <- 1
  }

  if(verbose) cat("Start_blockinfo_calculation\n")
  blockinfo <- list()
  indi <- ncol(dhm)
  if(max_groups==0){
    ite <- nrow(window_sequence)
  } else{
    ite <- nrow(dhm)
    window_sequence <- matrix(0, nrow=ite, ncol=4)
    window_sequence[1,] <- c(1,2,1+merging_error,1)
  }


  last.block <- rep(1, indi)
  nblocks.old <- 1
  old.ac <- 1
  if(max_groups==0){
    for(index in 1:ite){

      blockinfo[[index]] <- list()
      if(ite<=25 || index%%round(ite/25.1)==0){
        if(verbose) cat(".")
      }
      factor <- factorSNPs(c_dhm, window_sequence[index,1],
                           window_sequence[index,2])
      vers <- dhm[window_sequence[index,1]: window_sequence[index,2],
                  attr(factor, "where.to.find"), drop=FALSE]

      nblocks <- ncol(vers)
      similarity <- colSumsEqualSNPs(vers)

      if(blockinfo_mode==1){
        # Similarity-Merging
        skip <- rep(FALSE, nblocks)
        nr.coding <- 1:nblocks
        for(index3 in 1:nblocks){
          if(skip[index3]==FALSE){
            options <- similarity[index3,]
            merges <- unique(c(0,(options>=(window_sequence[index,4])) * (1:nblocks)))[-1]
            if(length(merges)>1){
              options <- similarity[merges,]
              merges <- unique(c(0,t(options>=(window_sequence[index,4])) * (1:nblocks)))[-1]
            }
            skip[merges] <- TRUE
            nr.coding[merges] <- index3
          }
        }

        current.block <- factor
        for(switch in 1:nblocks){
          if(nr.coding[switch]!=switch){
            current.block[current.block==switch] <- nr.coding[switch]
          }
        }
        hblocks <- numeric(nblocks)
        for(count in 1:nblocks){
          hblocks[count] <- sum(current.block==count)
        }
      } else{
        current.block <- factor
        hblocks <- attr(factor, "counts")
        ordering <- sort(hblocks, decreasing=TRUE, index.return=TRUE)$ix

        current <- 1
        skip <- rep(FALSE, nblocks)
        nr.coding <- 1:nblocks
        for(index3 in ordering){
          new_p <- ordering[which(similarity[ordering[1:current], index3]>=window_sequence[index,4])[1]]
          if(new_p!=index3){
            skip[index3] <- TRUE
            nr.coding[index3] <- new_p
            current.block[current.block==index3] <- nr.coding[new_p]
            similarity[index3,] <- similarity[new_p,]
            similarity[,index3] <- similarity[,new_p]
            current.block[current.block==index3] <- new_p
          }
          current <- current+1


        }

        hblocks <- numeric(nblocks)
        for(count in 1:nblocks){
          hblocks[count] <- sum(current.block==count)
        }
      }


      transition <- matrix(0, ncol=nblocks.old, nrow=nblocks)
      for(new in 1:nblocks){
        for(old in 1:nblocks.old){
          transition[new,old] <- sum((current.block==new)*(last.block==old))
        }
      }
      ac <- unique(nr.coding)
      blockinfo[[index]][[1]] <- hblocks[ac]
      blockinfo[[index]][[2]] <- similarity[ac,ac, drop=FALSE]
      blockinfo[[index]][[3]] <- transition[ac,old.ac, drop=FALSE]

      blockinfo[[index]][[4]] <- list()
      for(listn in 1:length(ac)){
        blockinfo[[index]][[4]][[listn]] <- vers[,nr.coding==ac[listn],drop=FALSE]
      }
      blockinfo[[index]][[5]] <- list()
      for(listn in 1:length(ac)){
        blockinfo[[index]][[5]][[listn]] <- which(current.block==ac[listn])
      }

      last.block <- current.block
      nblocks.old <- nblocks
      old.ac <- ac

      pos <- numeric(length(vers))
    }
  } else{
    activ_end <- 2
    current_window <- 1
    while(activ_end <= nrow(dhm)){
      possiblock <- list()
      if(ite<=25 || current_window%%round(ite/25.1)==0){
        if(verbose) cat(".")
      }


      factor <- factorSNPs(c_dhm, window_sequence[current_window,1],
                           window_sequence[current_window,2])
      vers <- dhm[window_sequence[current_window,1]: window_sequence[current_window,2],
                  attr(factor, "where.to.find"), drop=FALSE]

      nblocks <- ncol(vers)
      similarity <- colSumsEqualSNPs(vers)

      if(blockinfo_mode==1){
        skip <- rep(FALSE, nblocks)
        nr.coding <- 1:nblocks
        for(index3 in 1:nblocks){
          if(skip[index3]==FALSE){
            options <- similarity[index3,]
            merges <- unique(c(0,(options>=(window_sequence[current_window,4])) * (1:nblocks)))[-1]
            if(length(merges)>1){
              options <- similarity[merges,]
              merges <- unique(c(0,t(options>=(window_sequence[current_window,4])) * (1:nblocks)))[-1]
            }
            skip[merges] <- TRUE
            nr.coding[merges] <- index3
          }
        }

        current.block <- factor
        for(switch in 1:nblocks){
          if(nr.coding[switch]!=switch){
            current.block[current.block==switch] <- nr.coding[switch]
          }
        }
        hblocks <- numeric(nblocks)
        for(count in 1:nblocks){
          hblocks[count] <- sum(current.block==count)
        }
      } else{

        current.block <- factor
        hblocks <- attr(factor, "counts")
        ordering <- sort(hblocks, decreasing=TRUE, index.return=TRUE)$ix

        current <- 1
        skip <- rep(FALSE, nblocks)
        nr.coding <- 1:nblocks
        for(index3 in ordering){
          new_p <- ordering[which(similarity[ordering[1:current], index3]>=window_sequence[current_window,4])[1]]
          if(new_p!=index3){
            skip[index3] <- TRUE
            nr.coding[index3] <- new_p
            current.block[current.block==index3] <- nr.coding[new_p]
            similarity[index3,] <- similarity[new_p,]
            similarity[,index3] <- similarity[,new_p]
            current.block[current.block==index3] <- new_p
          }
          current <- current+1


        }

        hblocks <- numeric(nblocks)
        for(count in 1:nblocks){
          hblocks[count] <- sum(current.block==count)
        }
      }


      transition <- matrix(0, ncol=nblocks.old, nrow=nblocks)
      for(new in 1:nblocks){
        for(old in 1:nblocks.old){
          transition[new,old] <- sum((current.block==new)*(last.block==old))
        }
      }
      ac <- unique(nr.coding)
      possiblock[[1]] <- hblocks[ac]
      possiblock[[2]] <- similarity[ac,ac, drop=FALSE]
      possiblock[[3]] <- transition[ac,old.ac, drop=FALSE]
      possiblock[[4]] <- list()
      for(listn in 1:length(ac)){
        possiblock[[4]][[listn]] <- vers[,nr.coding==ac[listn], drop=FALSE]
      }
      possiblock[[5]] <- list()
      for(listn in 1:length(ac)){
        possiblock[[5]][[listn]] <- which(current.block==ac[listn])
      }
      if(activ_end==(nrow(dhm))){
        blockinfo[[current_window]] <- possiblock
        window_sequence <- window_sequence[1:current_window,]
        activ_end <- activ_end +1
      } else if(sum(hblocks>0)>max_groups){
        blockinfo[[current_window]] <- possiblock.temp
        last.block <- last.block.temp
        nblocks.old <- nblocks.old.temp
        old.ac <- old.ac.temp
        window_sequence[current_window,] <- window_sequence[current_window,] + c(0,-1,-1,-1)
        current_window <- current_window + 1
        activ_end <- activ_end + floor(log(max_groups, base=2)) - 1
        window_sequence[current_window, ] <- c(window_sequence[current_window-1, 2] +1 , activ_end, 0,0)
        window_sequence[current_window, 3] <- diff(window_sequence[current_window, 1:2]) + 1
        window_sequence[current_window, 4] <- window_sequence[current_window, 3] - merging_error
      } else{
        possiblock.temp <- possiblock
        last.block.temp <- current.block
        nblocks.old.temp <- nblocks
        old.ac.temp <- ac
        activ_end <- activ_end + 1
        window_sequence[current_window, 2:4] <- window_sequence[current_window, 2:4] +1
      }

    }
  }


  for(index in 1:length(blockinfo)){
    blockinfo[[index]][[6]] <- list()
    for(index2 in 1:length(blockinfo[[index]][[4]])){
      count <- numeric(ncol(blockinfo[[index]][[4]][[index2]]))
      if(length(count)>1){
        for(index3 in 1:length(count)){
          count[index3] <- sum(dhm[window_sequence[index,1]:window_sequence[index,2], blockinfo[[index]][[5]][[index2]]] == blockinfo[[index]][[4]][[index2]][,index3])
        }
      }else{
        count <- 1
      }
      blockinfo[[index]][[6]][[index2]] <- blockinfo[[index]][[4]][[index2]][,which.max(count)]

    }
  }
  if(length(bp_map)==0){
    bp_map <- rep(0, max(window_sequence[,2]))
  }
  for(index in 1:length(blockinfo)){
    blockinfo[[index]][[7]] <- c(bp_map[window_sequence[index,1]], bp_map[window_sequence[index,2]])
  }
  if(length(bp_map)>0 && ncol(window_sequence)==4){
    window_sequence <- cbind(window_sequence, bp_map[window_sequence[,1]], bp_map[window_sequence[,2]])
  }

  return(list(blockinfo, window_sequence))
}


