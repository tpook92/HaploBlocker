#' Extend-Block
#'
#' Function to add additional windows to each block
#' @param blocklist blocklist created by blocklist_calculcation algorithm
#' @param window_sequence_list sequence of predefined windows (default: NULL ;per row: start$snp, end$snp, length, length - merging_error, start$bp, end$bp)
#' @param dataset dataset which variant nr. for each window
#' @param indi number of haplotypes in the dataset
#' @param nwindow number of windows in the dataset
#' @param max_extending_diff Maximum number of windows with different realisation in the block-extending-algorithm
#' @param extending_ratio Minimum Ratio between windows with one different realisation to multiple in block-extending-algorithm
#' @return haplotype library
#'
extend_block <- function(blocklist, indi, nwindow, max_extending_diff=1, extending_ratio=10, dataset, window_sequence_list){
  prev <- numeric(max(nwindow))
  extensions_done <- 0
  for(index in 1:length(blocklist)){
    activ <- blocklist[[index]]

    #Verlaengerung nach Vorne
    if(activ[[2]]$window > 1 ){

      diff <- 0
      same <- 0
      extension_length <- 0
      max_diff <- 0
      position_diff <- NULL
      position <- activ[[2]]$window -1
      while(diff <= max_extending_diff && position>0){
        uniques <- unique(dataset[[activ[[12]]]][position, activ[[6]]])
        count <- numeric(length(uniques))
        for(index2 in 1:length(count)){
          count[index2] <- sum(dataset[[activ[[12]]]][position, activ[[6]]] == uniques[index2])
        }
        prev[position] <- uniques[which.max(count)]
        if(length(uniques)==1){
          same <- same +1
        } else{
          diff <- diff + 1
          position_diff <- c(position, position_diff)
        }
        if(same/diff >= max_diff && length(uniques) == 1){
          changes <- position_diff
          extension_length <- activ[[2]]$window - position
          max_diff <- same/diff
        }
        position <- position -1
      }

      if(max_diff > extending_ratio){
        extensions_done <- extensions_done +1
        activ[[2]]$window <- activ[[2]]$window - extension_length
        activ[[2]]$snp <- window_sequence_list[[activ[[12]]]][activ[[2]]$window,1]
        activ[[2]]$bp <- window_sequence_list[[activ[[12]]]][activ[[2]]$window,5]
        text <- paste("Extension_", activ[[2]]$window,"_", activ[[2]]$window + extension_length, sep="")
        activ[[1]] <- c(text, activ[[1]])
        activ[[4]] <- c(prev[activ[[2]]$window + 1:extension_length -1],activ[[4]])

        if(length(activ)>=8 && length(c(activ[[8]], changes))>0){
          activ[[8]] <- sort(c(activ[[8]], changes))
        } else if(length(changes)>0){
          activ[[8]] <- sort(changes)
        }
        blocklist[[index]] <- activ
      }

    }

    # Verlaengerung nach Hinten
    activ <- blocklist[[index]]
    if(activ[[3]]$window < nwindow[[activ[[12]]]] ){

      diff <- 0
      same <- 0
      extension_length <- 0
      max_diff <- 0
      position_diff <- NULL
      position <- activ[[3]]$window +1
      while(diff <= max_extending_diff && position<= nwindow[[activ[[12]]]]){
        uniques <- unique(dataset[[activ[[12]]]][position, activ[[6]]])
        count <- numeric(length(uniques))
        for(index2 in 1:length(count)){
          count[index2] <- sum(dataset[[activ[[12]]]][position, activ[[6]]] == uniques[index2])
        }
        prev[position] <- uniques[which.max(count)]
        if(length(uniques)==1){
          same <- same +1
        } else{
          diff <- diff + 1
          position_diff <- c(position_diff, position)
        }
        if(same/diff >= max_diff && length(uniques) == 1){
          changes <- position_diff
          extension_length <- position - activ[[3]]$window
          max_diff <- same/diff
        }
        position <- position +1
      }

      if(max_diff > extending_ratio){
        extensions_done <- extensions_done +1
        activ[[3]]$window <- activ[[3]]$window + extension_length
        activ[[3]]$snp <- window_sequence_list[[activ[[12]]]][activ[[3]]$window,2]
        activ[[3]]$bp <- window_sequence_list[[activ[[12]]]][activ[[3]]$window,6]
        text <- paste("Extension_", activ[[3]]$window- extension_length,"_", activ[[3]]$window , sep="")
        activ[[1]] <- c(activ[[1]], text)
        activ[[4]] <- c(activ[[4]],  prev[activ[[3]]$window + 1:extension_length - extension_length])
        if(length(activ)>=8 && length(c(activ[[8]], changes))>0){
          activ[[8]] <- sort(c(activ[[8]], changes))
        } else if(length(changes)>0){
          activ[[8]] <- sort(changes)
        }
        blocklist[[index]] <- activ

      }

    }



  }

  return(list(blocklist, extensions_done))
}
