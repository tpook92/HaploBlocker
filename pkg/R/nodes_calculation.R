#' Calculate transition probabilities between nodes
#'
#' Internal Function to calculate calculate a node-dataset
#' @param blockinfo List with all relevant information to each window seperatly
#' @param window_sequence sequence of predefined windows (default: NULL ;per row: start$snp, end$snp, length, length - merging_error, start$bp, end$bp)
#' @param verbose Set to FALSE to not display any prints
#' @return window cluster without transition probabilities

nodes_calculation <- function(blockinfo, window_sequence, verbose = TRUE){
  if(verbose) cat("Start_nodes_calculation\n")
  nanimals <- sum(blockinfo[[1]][[1]])
  data <- list()
  nr <- 1
  n <- length(blockinfo)
  for(index in 1:n){
    if(n<=25 || index%%round(n/25.1)==0){
      cat(".")
    }
    for(index2 in 1:length(blockinfo[[index]][[1]])){
      start_l <- list()
      start_l$window <- index
      start_l$snp <- window_sequence[index, 1]
      start_l$bp <- blockinfo[[index]][[7]][1]
      end_l <- list()
      end_l$window <- index
      end_l$snp <- window_sequence[index,2]
      end_l$bp <- blockinfo[[index]][[7]][2]
      data[[nr]] <- list()
      data[[nr]]$start <- start_l
      data[[nr]]$end <- end_l
      data[[nr]]$n_haplo <- blockinfo[[index]][[1]][index2]
      data[[nr]]$window_sequence <- index2
      data[[nr]]$haplos <- blockinfo[[index]][[5]][[index2]]
      nr <- nr + 1

    }

  }
  return(data)
}




