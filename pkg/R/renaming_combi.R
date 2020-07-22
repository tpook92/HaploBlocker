#' Renaming function
#'
#' Internal Function to remove empty nodes and rename after merging
#' @param data1 node-dataset
#' @param nwindow number of windows in the dataset
#' @return Updated window cluster

renaming_combi <- function(data1, nwindow){
  b <- nodes_size(data1)
  a <- start_end_block(data1)

  relevant <- which(b>0)
  if(length(relevant)==0){
    print("Unable to identify common variants in the dataset - please reduce window_size or allow for less haplotypes per node (node_min)!")
  }
  pos <-sort(a[relevant,1], index.return=TRUE)$ix
  old.name <- relevant[pos]

  new.name.code <- numeric(max(old.name))
  new.name.code[old.name] <- 1:length(old.name)

  new.data <- list()
  new.length <- sum(b>0)


  nmax <- length(old.name)

  for(index in 1:new.length){
    new.data[[index]] <- data1[[old.name[index]]]
    #switchname
    if(new.data[[index]][[2]]$window < nwindow){
      old_a <- new.data[[index]][[6]][,1]
      row <- 1
      for(switch in old_a[old_a!=0]){
        new.data[[index]][[6]][row,1] <- new.name.code[switch]
        row <- row+1
      }
    }

    if(new.data[[index]][[1]]$window>1){
      old_e <- new.data[[index]][[7]][,1]
      row <- 1
      for(switch in old_e[old_e!=0]){

        new.data[[index]][[7]][row,1] <- new.name.code[switch]


        row <- row+1
      }
    }

  }
  return(new.data)
}
