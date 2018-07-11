#' Mainfunction to calculate haplotype blocks
#'
#' Function to generate haplotype blocks from haploid data
#' @param dhm haploid SNP-dataset
#' @param window_size size of each window in the algorithm (default: 20)
#' @param merging_error number of allowed errors per block (default: 1)
#' @param window_sequence sequence of predefined windows (default: NULL ;per row: start$snp, end$snp, length, length - merging_error, start$bp, end$bp)
#' @param bp_map vector of positions for each SNP in bp (default: NULL - all 0)
#' @param window_anchor_gens matrix to constructed window_sequence base on start/end points in bp (e.g. gen regions, per row: start, gen)
#' @param at_least_one If TRUE no allowed merging errors in windows of size 1
#' @param blockinfo_mode Structure of the groups in step I (default: 0 - Common haplos as major variants, 1- minimum number of groups)
#' @param actual_snp_weight Set weight for difference between two alleles in a SNP (more than 1 possible base pair)
#' @param na_snp_weight Set weight for difference between NA and allele in a SNP (more than 1 possible base pair)
#' @param na_seq_weight Set weight for difference between NA and allele a loci with 1 possible base pair
#' @export


blockinfo_calculation_na <- function(dhm, window_sequence=NULL, window_size=32, merging_error=3, bp_map=NULL,
                                  window_anchor_gens=NULL, at_least_one=TRUE, blockinfo_mode=0, actual_snp_weight = 5,
                                  na_snp_weight=2, na_seq_weight=0){

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

  if(length(window_sequence)==0){
    window_sequence <- cbind(0:(nrow(dhm)/window_size-1)*window_size+1, 1:(nrow(dhm)/window_size)*window_size)
  }
  if(ncol(window_sequence)==2){
    window_sequence <- cbind(window_sequence, window_sequence[,2]-window_sequence[,1]+1, window_sequence[,2]-window_sequence[,1]-merging_error+1)
  }
  if(ncol(window_sequence)==3){
    window_sequence <- cbind(window_sequence[,1:2], window_sequence[,2]-window_sequence[,1]+1, window_sequence[,2]-window_sequence[,1]+1-window_sequence[,3])
  }
  if(at_least_one== TRUE && sum(window_sequence[,4]==0)>0){
    window_sequence[window_sequence[,4]==0, 4] <- 1
  }

  cat("Start_blockinfo_calculation\n")
  blockinfo <- list()
  indi <- ncol(dhm)
  ite <- nrow(window_sequence)

  last.block <- rep(1, indi)
  nblocks.old <- 1
  old.ac <- 1
  for(index in 1:ite){
    blockinfo[[index]] <- list()
    if(ite<=25 || index%%round(ite/25.1)==0){
      cat(".")
    }
    activ <- dhm[window_sequence[index,1]:window_sequence[index,2],]
    if(is.matrix(activ)==FALSE){
      activ <- t(as.matrix(activ))
    }
    snp_weight <- rep(1, window_sequence[index,3])
    na_weight <- rep(na_seq_weight, window_sequence[index,3])
    if(TRUE){
      for(checkor in 1:window_sequence[index,3]){
        possible_snp <- unique(c(9, activ[checkor,]))
        if(length(possible_snp)>2){
          snp_weight[checkor] <- actual_snp_weight
          na_weight[checkor] <- actual_snp_weight - na_snp_weight
          window_sequence[index,4] <- window_sequence[index,4] +actual_snp_weight -1
        }
      }
    }

    vers <- t(unique(activ, MARGIN=2))
    nblocks <- nrow(vers)

    similarity <- matrix(0, ncol=nblocks, nrow=nblocks)
    for(col in 1:nblocks){
      for(row in 1:col){

        similarity[row,col] <- sum((vers[col,]==vers[row,])*snp_weight) + sum( ((vers[col,]==9)*(vers[row,]!=9) + (vers[col,]==9)*(vers[row,]!=9)) * na_weight)
        similarity[col,row] <- similarity[row,col]
      }
    }

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


      options <- numeric(nrow(vers))
      for(index2 in 1:length(options)){
        options[index2] <- paste0(vers[index2,], collapse = "")
      }
      func1 <- function(x){
        x0 <- paste(x, collapse="")
        which(options==x0)
      }
      current.block <- apply(activ, MARGIN=2, func1)


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


      options <- numeric(nrow(vers))
      for(index2 in 1:length(options)){
        options[index2] <- paste0(vers[index2,], collapse = "")
      }
      func1 <- function(x){
        x0 <- paste(x, collapse="")
        which(options==x0)
      }
      current.block <- apply(activ, MARGIN=2, func1)


      hblocks <- numeric(nblocks)
      for(count in 1:nblocks){
        hblocks[count] <- sum(current.block==count)
      }

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
    blockinfo[[index]][[2]] <- rbind(similarity[ac,ac],NULL)
    blockinfo[[index]][[3]] <- transition[ac,old.ac]
    if(is.matrix(transition[ac,old.ac])==FALSE){
      if(length(old.ac)==1) {
        blockinfo[[index]][[3]] <- cbind(transition[ac,old.ac], NULL)
      } else{
        blockinfo[[index]][[3]] <- matrix(transition[ac,old.ac], nrow=1)
      }
    }
    blockinfo[[index]][[4]] <- list()
    for(listn in 1:length(ac)){
      blockinfo[[index]][[4]][[listn]] <- rbind(vers[nr.coding==ac[listn],],NULL)
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

  for(index in 1:length(blockinfo)){
    blockinfo[[index]][[6]] <- list()
    for(index2 in 1:length(blockinfo[[index]][[4]])){
      count <- numeric(nrow(blockinfo[[index]][[4]][[index2]]))
      if(length(count)>1){
        for(index3 in 1:length(count)){
          count[index3] <- sum(dhm[window_sequence[index,1]:window_sequence[index,2], blockinfo[[index]][[5]][[index2]]] == blockinfo[[index]][[4]][[index2]][index3,])
        }
      }else{
        count <- 1
      }
      blockinfo[[index]][[6]][[index2]] <- blockinfo[[index]][[4]][[index2]][which.max(count),]

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


