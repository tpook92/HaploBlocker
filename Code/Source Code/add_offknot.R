#' Off-variant-identification
#'
#' Function to add off-variant-blocks
#' @param blocklist blocklist created by blocklist_calculcation algorithm
#' @param window_sequence sequence of predefined windows (default: NULL ;per row: start$snp, end$snp, length, length - merging_error, start$bp, end$bp)
#' @param dataset dataset which variant nr. for each window
#' @param indi number of haplotypes in the dataset
#' @param nwindow number of windows in the dataset
#' @param bp_map vector of positions for each SNP in bp (default: NULL - all 0)
#' @param off_knot_minimum_blocklength Minimum length of newly identified blocks (default: 10)
#' @param off_knot_minimum_blocksize Minimum number of individuals in newly identified blocks (default: 5)
#' @param raster Raster-width in the identification step (default: 5; recommended to be lower than off_knot_minimum_blocklength)
#' @export

add_offknot <- function(blocklist, dataset, indi, nwindow, window_sequence, bp_map, off_knot_minimum_blocklength=10, off_knot_minimum_blocksize=5, raster=5){
  t <- coverage_test(blocklist, indi, type="window", max=1)
  index <- 1
  while(index+off_knot_minimum_blocklength -1 <= nwindow){
    window <- index:(index+off_knot_minimum_blocklength-1)
    center <- round(index+off_knot_minimum_blocklength/2)
    non_block_indi <- unique(c(0,(t[center,]==0) * (1:ncol(t))))[-1]
    activ_data <- dataset[window, non_block_indi]
    varants <- unique(t(activ_data))
    counts <- numeric(nrow(varants))
    for(index2 in 1:length(counts)){
      counts[index2] <- sum(colSums(activ_data==varants[index2,])==off_knot_minimum_blocklength)
    }
    varants <- varants[counts>=off_knot_minimum_blocksize,,drop=FALSE]

    if(nrow(varants)>0){
      for(index2 in 1:nrow(varants)){
        activ_indi <- non_block_indi[which(colSums(activ_data==varants[index2,])==off_knot_minimum_blocklength)]
        start <- index+1
        end <- index+off_knot_minimum_blocklength-1 -1
        extents <- TRUE
        extentb <- TRUE
        while(extents==TRUE && start>1){
          start <- start - 1
          if(length(unique(dataset[start-1, activ_indi]))!=1){
            extents <- FALSE
          }
        }
        while(extentb==TRUE && end<(nwindow-1)){
          end <- end +1
          if(length(unique(dataset[end+1, activ_indi]))!=1){
            extentb <- FALSE
          }
        }
        end <- end + extentb # NUR TRUE wenn am ende des Datensatzes

        new_block <- list()

        new_block[[1]] <- "Off-Window-Variant"
        new_block[[2]] <- list()
        new_block[[2]]$window <- start
        new_block[[2]]$snp <- window_sequence[start,1]

        new_block[[3]] <- list()
        new_block[[3]]$window <- end
        new_block[[3]]$snp <- window_sequence[end,2]
        if(length(bp_map)==0){
          new_block[[2]]$bp <- 0
          new_block[[3]]$bp <- 0
        } else{
          new_block[[2]]$bp <- bp_map[new_block[[2]]$snp]
          new_block[[3]]$bp <- bp_map[new_block[[3]]$snp]
        }

        new_block[[4]] <- dataset[start:end, activ_indi[1]]
        new_block[[5]] <- length(activ_indi)
        new_block[[6]] <- activ_indi

        blocklist[[length(blocklist)+1]] <- new_block
      }
    }




    index <- index+raster
  }
  return(blocklist)
}
