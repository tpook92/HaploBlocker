#' Renaming function
#'
#' Internal Function to remove empty nodes and rename after merging
#' @param data node-dataset
#' @param nwindow number of windows in the dataset
#' @export

renaming_combi <- function(data, nwindow){
  b <- nodes_size(data)
  a <- start_end_block(data)

  relevant <- which(b>0)
  pos <-sort(a[relevant,1], index.return=TRUE)$ix
  old.name <- relevant[pos]

  new.data <- list()
  new.length <- sum(b>0)
  for(index in 1:new.length){
    new.data[[index]] <- data[[old.name[index]]]
    #switchname
    if(new.data[[index]][[2]]$window < nwindow){
      old_a <- new.data[[index]][[6]][,1]
      row <- 1
      for(switch in old_a[old_a!=0]){

        new.data[[index]][[6]][row,1] <- which(old.name==switch)
        row <- row+1
      }
    }

    if(new.data[[index]][[1]]$window>1){
      old_e <- new.data[[index]][[7]][,1]
      row <- 1
      for(switch in old_e[old_e!=0]){

        new.data[[index]][[7]][row,1] <- which(old.name==switch)
        row <- row+1
      }
    }


  }
  return(new.data)
}
