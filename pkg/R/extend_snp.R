#' Extend-SNP
#'
#' Function to add single SNPs to each block
#' @param dhm haploid SNP-dataset
#' @param blocklist blocklist created by blocklist_calculcation algorithm
#' @param window_sequence_list sequence of predefined windows (default: NULL ;per row: start$snp, end$snp, length, length - merging_error, start$bp, end$bp)
#' @param indi number of haplotypes in the dataset
#' @param nwindow number of windows in the dataset
#' @param bp_map vector of positions for each SNP in bp (default: NULL - all 0)
#' @param max_extending_diff_snp Maximum number of SNPs with variants in SNP-extending-algorithm (step V; default: 0)
#' @param extending_ratio_snp Minimum ratio of SNPs with only one allele to those with variants (default: Inf)#' @param off_node_addition If TRUE use off-variant-identification (default: FALSE)
#' @export


extend_snp <- function(blocklist, indi, nwindow, dhm, window_sequence_list, bp_map, max_extending_diff_snp=0, extending_ratio_snp=Inf){
  max_l <- 0
  for(index in 1:length(window_sequence_list)){
    if(is.matrix(window_sequence_list[[index]])){
      max_l <- max(max_l, window_sequence_list[[index]][,3])
    }

  }
  prev <- numeric(max_l)
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
      position <- activ[[2]]$snp -1
      while(diff <= max_extending_diff_snp && position>0){
        count <- sum(dhm[position, activ[[6]]]==dhm[position, 1])
        if(count==activ[[5]] || count==0){
          same <- same +1
          prev[same+diff] <- dhm[position, activ[[6]][1]]
        } else{
          diff <- diff + 1
          position_diff <- c(position, position_diff)
        }
        if(same/diff >= max_diff){
          extension_length <- activ[[2]]$snp - position
          max_diff <- same/diff
        }
        position <- position -1
      }

      if(max_diff >= extending_ratio_snp){
        if(extension_length==max_l){
          cat("SNP-extension over full window. That is not supposed to happen\n")
        }
        extensions_done <- extensions_done +1
        activ[[2]]$snp <- activ[[2]]$snp - extension_length
        activ[[2]]$bp <- bp_map[activ[[2]]$snp]
        text <- paste("SNP-Extension:",activ[[2]]$snp,"_",  activ[[2]]$snp + extension_length-1, sep="")
        activ[[1]] <- c(text, activ[[1]])
        #activ[[4]] <- c(prev[activ[[2]]$window + 1:extension_length -1],activ[[4]])

        blocklist[[index]] <- activ
      }

    }

    # Verlaengerung nach Hinten
    activ <- blocklist[[index]]
    if(activ[[3]]$window < nwindow[activ[[12]]] ){

      diff <- 0
      same <- 0
      extension_length <- 0
      max_diff <- 0
      position_diff <- NULL
      position <- activ[[3]]$snp +1
      while(diff <= max_extending_diff_snp && position<= max(window_sequence_list[[activ[[12]]]][,2])){
        count <- sum(dhm[position, activ[[6]]]==dhm[position, 1])
        if(count==activ[[5]] || count==0){
          same <- same +1
          prev[same+diff] <- dhm[position, activ[[6]][1]]
        } else{
          diff <- diff + 1
          position_diff <- c(position, position_diff)
        }
        if(same/diff >= max_diff){
          extension_length <- position - activ[[3]]$snp
          max_diff <- same/diff
        }

        position <- position +1
      }

      if(max_diff >= extending_ratio_snp){
        if(extension_length==max_l){
          cat("SNP-extension over full window. That is not supposed to happen\n")
        }
        extensions_done <- extensions_done +1
        activ[[3]]$snp <- activ[[3]]$snp + extension_length
        activ[[3]]$bp <- bp_map[activ[[3]]$snp]
        text <- paste("SNP-Extension:",activ[[2]]$snp - extension_length+1,"_",  activ[[2]]$snp, sep="")
        activ[[1]] <- c(activ[[1]], text)
        #activ[[4]] <- c(activ[[4]],  prev[activ[[3]]$window + 1:extension_length - extension_length])
        blocklist[[index]] <- activ

      }

    }

  }

  return(blocklist)
}
