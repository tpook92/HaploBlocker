#' Determine number of differnt haplotypes per window
#'
#' Function to determine number of differnt haplotypes per window
#' @param blockinfo List with all relevant information to each window seperatly
#' @return Number of variants in each window

blockinfo_size <- function(blockinfo){
  size <- numeric(length(blockinfo))
  for(index in 1:length(blockinfo)){
    size[index] <- length(blockinfo[[index]][[1]])
  }
  return(size)
}

#' Determine number of differnt haplotypes per window
#'
#' Function to determine number of differnt haplotypes per window (before error reduction)
#' @param blockinfo List with all relevant information to each window seperatly
#' @return Number of variants (before merging) in each window

blockinfo_size0 <- function(blockinfo){
  size <- numeric(length(blockinfo))
  for(index in 1:length(blockinfo)){
    for(index2 in 1:length(blockinfo[[index]][[1]])){
      size[index] <- size[index] + nrow(blockinfo[[index]][[4]][[index2]])
    }
  }
  return(size)
}

#' Determine the maximum number of the same haplotype in each window
#'
#' Function to determine the number of the same haplotype in each window
#' @param blockinfo List with all relevant information to each window seperatly
#' @return Maximum number of same variants in each window

blockinfo_max <- function(blockinfo){
  size <- numeric(length(blockinfo))
  for(index in 1:length(blockinfo)){
    size[index] <- max(blockinfo[[index]][[1]])
  }
  return(size)
}
