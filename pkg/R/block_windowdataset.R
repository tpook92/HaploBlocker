#' Function to derive block-dataset (uniform-window)
#'
#' Function to compute a block-dataset with uniform start/end points
#' @param blocklist blocklist
#' @param data haploid SNP-dataset
#' @param consider_nonblock Consider haplotypes that are in no block (default: FALSE)
#' @param return_dataset Return a per window dataset with more than 2 variants instead of binary coded one (default: FALSE)
#' @param non_haploblocker Set TRUE to consider all haplotypes and all variants in each window
#' @param non_haploblocker_merging_error Set the number of markers with potential variations from the main variant (default: 0)
#' @param min_freq_boundary Set the minimum frequency of a block in the population (default: 0)
#' @param min_length_window Set the minimum length of each window (default: 1)
#' @param verbose Set to FALSE to not display any prints
#' @examples
#' data(blocklist_ex_maze)
#' windowdataset <- block_windowdataset(blocklist_ex_maze)
#' @export
#' @return Windowblock based dataset


block_windowdataset <- function(blocklist=NULL, data=NULL, consider_nonblock=FALSE, return_dataset=FALSE, non_haploblocker=FALSE,
                                non_haploblocker_merging_error=0, min_freq_boundary=0, min_length_window=1, verbose= TRUE){

  if(length(data)>0 && is.matrix(data)==FALSE){
    data <- as.matrix(data)
  }
  if(length(blocklist)==0){
    blocklist <- block_calculation(data, verbose= verbose)
  }
  n_indi <- indi_calc(blocklist)



  se <- blocklist_startend(blocklist)

  end_block  <- sort(unique(c(0,se[,1]-1, se[,2])))[-1]

  if(min_freq_boundary>0){
    if(min_freq_boundary<1){
      min_freq_boundary <- n_indi * min_freq_boundary
    }

    ## Freq
    freq_end <- numeric(length(end_block))
    vars <- list()
    for(index in 1:length(blocklist)){
      count <- which((blocklist[[index]][[2]]$snp-1) == end_block)
      if(length(count)>0){
        if(length(vars)>=count){
          vars[[count]] <- unique(c(vars[[count]], blocklist[[index]][[6]]))
        } else{
          vars[[count]] <- blocklist[[index]][[6]]
        }
      }

      count <- which((blocklist[[index]][[3]]$snp) == end_block)
      if(length(count)>0){
        if(length(vars)>=count){
          vars[[count]] <- unique(c(vars[[count]], blocklist[[index]][[6]]))
        } else{
          vars[[count]] <- blocklist[[index]][[6]]
        }
      }


    }
    for(index in 1:length(vars)){
      freq_end[index] <- length(vars[[index]])
    }

    remove <- which(freq_end < min_freq_boundary)
    end_block <- end_block[-remove]

  }

  start_block <- c(min(se), end_block[1:(length(end_block)-1)]+1)

  if(min_length_window>0){

    block_length <- end_block - start_block + 1
    merging_candidate <- which(block_length<min_length_window)
    while(length(merging_candidate)>0){
      for(index in merging_candidate){
        if(min(block_length[index + c(-2,-1,0,1,2)])==block_length[index]){

          if(index==length(block_length) || (index >1 && (block_length[index-1]<block_length[index+1]))){
            end_block[index-1] <- end_block[index]
            block_length[index-1] <- end_block[index-1] - start_block[index-1] +1
          } else{
            start_block[index+1] <- start_block[index]
            block_length[index+1] <- end_block[index+1] - start_block[index+1] +1
          }
          start_block[index] <- 0
          end_block[index] <- 0
          block_length[index] <- -Inf

        }

      }
      remove <- which(start_block==0)
      end_block <- end_block[-remove]
      start_block <- start_block[-remove]
      block_length <- block_length[-remove]

      merging_candidate <- which(block_length<min_length_window)
    }


  }
  if(non_haploblocker){
    dhm <- as.matrix(data)
    unique.dhm <- unique(as.vector(dhm))
    fixcoding(unique.dhm)
    c_dhm <- codeSNPs(dhm)
    blockinfo <- blockinfo_calculation(dhm, window_sequence = cbind(start_block, end_block), merging_error = non_haploblocker_merging_error, c_dhm=c_dhm,
                                       verbose=verbose)
    dataset <- block_dataset_construction(blockinfo[[1]],blocklist)
  } else{
    dataset <- matrix(0, nrow=length(start_block), ncol=n_indi)

    for(index in 1:length(start_block)){
      active_blocks <- which(((se[,1]<=start_block[index]) * (se[,2]>=end_block[index])==1))
      if(length(active_blocks)>0){
        window <- c(start_block[index], end_block[index])
        variants <- matrix(NA, nrow=length(active_blocks), ncol=diff(window)+1)
        for(index2 in 1:length(active_blocks)){
          variants[index2,] <- blocklist[[active_blocks[[index2]]]][[7]]$snp[-blocklist[[active_blocks[[index2]]]][[2]]$snp +window[1]:window[2]+1]
        }
        variants <- unique(variants)
        class <- list()
        class_count <- numeric(nrow(variants))

        for(index2 in 1:length(active_blocks)){
          class_nr <- which.max(colSums(t(variants)==blocklist[[active_blocks[[index2]]]][[7]]$snp[-blocklist[[active_blocks[[index2]]]][[2]]$snp +window[1]:window[2]+1]))
          if(length(class)<class_nr){
            class[[class_nr]] <- blocklist[[active_blocks[[index2]]]][[6]]
          } else{
            class[[class_nr]] <- unique(c(class[[class_nr]],blocklist[[active_blocks[[index2]]]][[6]]))
          }

        }

        for(index2 in 1:nrow(variants)){
          if(length(class)>=index2){
            class_count[index2] <- length( class[[index2]])
          }

        }

        while(sum(duplicated(unlist(class))>0)){

          full <- unlist(class)
          full <- full[which(duplicated(full))]
          candidates_remove <- NULL
          for(index2 in 1:length(class)){
            if(sum(class[[index2]]==full[1])){
              candidates_remove <- c(candidates_remove, index2)
            }
          }
          candidates_remove <- candidates_remove[-which.max(class_count[candidates_remove])]
          for(index2 in candidates_remove){
            class[[index2]] <- setdiff(class[[index2]], full)
          }
          for(index2 in 1:nrow(variants)){
            if(length(class)>=index2){
              class_count[index2] <- length( class[[index2]])
            }

          }

        }

        running <- 1
        class_index <- numeric(length(class))
        for(index2 in 1:length(class)){
          dataset[index, class[[index2]]] <- running
          if(length(class[[index2]])>0){
            class_index[index2] <- running
            running <- running + 1
          }
        }

        if(consider_nonblock){
          to_test <- which(dataset[index,]==0)
          for(index2 in to_test){
            temp1 <- colSums(data[window[1]:window[2], index2] == t(variants))
            if(max(temp1)==ncol(variants)){
              dataset[index,index2] <- class_index[which.max(temp1)]
            }
          }
        }
      }

    }

  }
  var_permarker <- RandomFieldsUtils::colMax(t(dataset))
  bin_dataset <- matrix(0, nrow=sum(var_permarker), ncol=n_indi)
  counter <- 1
  row_names <- numeric(nrow(bin_dataset))
  for(index in 1:length(var_permarker)){
    if(length(var_permarker[index])>0){
      if(var_permarker[index]>0){
        for(index2 in 1:var_permarker[index]){
          bin_dataset[counter,] <- dataset[index,]==index2
          row_names[counter] <- paste0("window:",start_block[index],"-", end_block[index], "variant", index2)
          counter <- counter + 1
        }
      }

    }
  }
  if(return_dataset){
    rownames(dataset) <- paste0("window:",start_block,"-", end_block)
    colnames(dataset) <- paste0("haplo", 1:ncol(dataset))
    return(dataset)
  } else{
    rownames(bin_dataset) <- row_names
    colnames(bin_dataset) <- paste0("haplo", 1:ncol(dataset))
    return(bin_dataset)
  }

}



