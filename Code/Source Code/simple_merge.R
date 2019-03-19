#' Fast Simple Merge nodes
#'
#' Faster version of simple_merge_probv2 to perform simple merge for node data
#' @param data node-dataset
#' @param intersect_func Used intersect-function (internally relevant for computation time)
#' @export


simple_merge <- function(data, intersect_func=intersect){
  cat(paste0("Start_simple_merge: ", length(data),"\n"))
  a <- start_end_block(data)
  relevant <- (a[,1]>1)
  relevant <- (1:nrow(a))[relevant]
  relevant <- relevant[length(relevant):1]
  prev <- Inf
  size <- data_size(data)
  for(check in relevant){
    activ <- data[[check]]
    if(prev!=activ[[1]]$window){
      possible <- which(a[,2]==(activ[[1]]$window-1))
      prev <- activ[[1]]$window
    }

    possible2 <- possible[size[possible]==size[check]]

    if(length(possible2)>0){
      for(index in 1:length(possible2)){
        transi <- length(intersect_func(data[[possible2[index]]][[5]], activ[[5]]))
        if(transi == data[[possible2[index]]][[3]]){
          data[[possible2[index]]][[2]] <- activ[[2]]
          data[[possible2[index]]][[4]] <- c(data[[possible2[index]]][[4]],activ[[4]])
          data[[check]] <- NULL
          break
        }
      }
    }

  }
  return(data)
}

