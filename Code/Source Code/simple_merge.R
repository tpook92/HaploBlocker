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
  for(check in relevant){

    activ <- data[[check]]
    possible <- (1:nrow(a))[(a[,2]==activ[[1]]$window-1)]
    for(index in 1:length(possible)){
      transi <- length(intersect_func(data[[possible[index]]][[5]], activ[[5]]))
      if(transi== length(data[[possible[index]]][[5]]) && transi== length(activ[[5]])){
        data[[possible[index]]][[2]] <- activ[[2]]
        data[[possible[index]]][[4]] <- c(data[[possible[index]]][[4]],activ[[4]])
        data[[check]] <- NULL
      }
    }
  }
  return(data)
}
