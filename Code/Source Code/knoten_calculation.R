#' Calculate transition probabilities between knots
#'
#' Internal Function to calculate calculate a knot-dataset
#' @param blockinfo List with all relevant information to each window seperatly
#' @param window_sequence sequence of predefined windows (default: NULL ;per row: start$snp, end$snp, length, length - merging_error, start$bp, end$bp)
#' @export

knoten_calculation <- function(blockinfo, window_sequence){
  print("Start_knoten_calculation")
  nanimals <- sum(blockinfo[[1]][[1]])
  next.block <- rep(1,nanimals)
  data <- list()
  max.block <- numeric(10000)
  nr <- 1
  n <- length(blockinfo)
  news <- numeric(n)
  for(index in 1:n){
    if(n<=25 || index%%round(n/25.1)==0){
      cat(".")
    }
    col <- ncol(blockinfo[[index]][[3]])
    maxi <- numeric(col)
    for(index2 in 1:col){
      maxi[index2] <- max(blockinfo[[index]][[3]][,index2])
    }
    if(index>1){
      news[index] <- length(blockinfo[[index-1]][[1]]) - sum(maxi == colSums(blockinfo[[index]][[3]]))
    } else{
      news[index] <- length(blockinfo[[1]][[1]])
    }


    if(news[index]>0){
      if(index==1){
        relevant <- 1
      } else{
        relevant <- unique(c(0,1:length(maxi) * (1-(maxi == colSums(blockinfo[[index]][[3]])))))[-1]
      }
      for(abc in relevant){
        varianten <- unique(c(0,(blockinfo[[index]][[3]][,abc]>0) * (1:nrow(blockinfo[[index]][[3]]))))[-1]
        for(dfg in varianten){
          start <- index
          anzahl <- blockinfo[[index]][[3]][dfg,abc]
          max.block[1] <- dfg
          blocklength <- 1
          while(index+blocklength <= n && sum(blockinfo[[index+blocklength]][[3]][,max.block[blocklength]])==max(blockinfo[[index+blocklength]][[3]][,max.block[blocklength]]) ){

            max.block[blocklength+1] <- which.max(blockinfo[[index+blocklength]][[3]][,max.block[blocklength]])
            blocklength <- blocklength + 1

          }
          end <- start + blocklength -1
          one <- numeric(nanimals)
          two <- numeric(nanimals)
          one[blockinfo[[start]][[5]][[dfg]]] <- 1
          if(index==1){
            two <- rep(1, nanimals)
          } else{
            two[blockinfo[[start-1]][[5]][[abc]]] <- 1
          }
          activ.animals <- unique(c(0,(1:nanimals) * (next.block==start) * one * two))[-1]
          next.block[activ.animals] <- end +1
          start_l <- list()
          start_l$window <- start
          start_l$snp <- window_sequence[start, 1]
          start_l$bp <- blockinfo[[start]][[7]][1]
          end_l <- list()
          end_l$window <- end
          end_l$snp <- window_sequence[end,2]
          end_l$bp <- blockinfo[[end]][[7]][2]
          data[[nr]] <- list()
          data[[nr]]$start <- start_l
          data[[nr]]$end <- end_l
          data[[nr]]$n_haplo <- anzahl
          data[[nr]]$window_sequence <- max.block[1:blocklength]
          data[[nr]]$haplos <- activ.animals
          nr <- nr + 1
        }
      }
    }
  }
  return(data)
}




