#' Calculate transition probabilities between nodes
#'
#' Internal Function to calculate transition probabilities between nodes
#' @param data1 node-dataset
#' @param intersect_func Used intersect-function (internally relevant for computation time)
#' @return window cluster


calculate_transition <- function(data1, intersect_func=intersect){

  a <- start_end_block(data1)
  relevant <- (a[,1]>1)
  relevant <- (1:nrow(a))[relevant]
  relevant <- relevant[length(relevant):1]
  for(index in which(a[,2]==max(a[,2]))){
    data1[[index]][[6]] <- cbind(0,0)
  }

  lastw <- 0
  ttt <- 1:nrow(a)

  for(check in relevant){
    activ <- data1[[check]]
    if(lastw != activ[[1]]$window){
      possible <- ttt[(a[,2]==activ[[1]]$window-1)]
      lastw <- activ[[1]]$window
    }


    for(index in 1:length(possible)){
      transi <- length(intersect(data1[[possible[index]]][[5]], activ[[5]]))
      if(transi>0){
        if(length(data1[[possible[index]]])==5){
          data1[[possible[index]]][[6]] <- cbind(check,transi)
        } else{
          data1[[possible[index]]][[6]] <- rbind(data1[[possible[index]]][[6]], c(check,transi) )
        }

      }


    }
  }


  ## Eingaenge

  relevant <- (a[,2]<max(a[,2]))
  relevant <- (1:nrow(a))[relevant]

  for(index in which(a[,1]==1)){
    data1[[index]][[7]] <- cbind(0,0)
  }

  for(check in relevant){
    nr <- 1
    poss <- data1[[check]][[6]][,1]
    for(index in poss){
      if(length(data1[[index]])<7){
        data1[[index]][[7]] <- matrix(c(check, data1[[check]][[6]][nr,2]), nrow=1)
      } else{
        data1[[index]][[7]] <- rbind( data1[[index]][[7]], c(check, data1[[check]][[6]][nr,2]))
      }

      nr <- nr +1
    }
  }
  return(data1)

}

'

calculate_transition <- function(data){

  a <- start_end_block(data)
  relevant <- (a[,1]>1)
  relevant <- (1:nrow(a))[relevant]
  relevant <- relevant[length(relevant):1]
  for(check in relevant){
    activ <- data[[check]]
    possible <- (1:nrow(a))[(a[,2]==activ[[1]]-1)]

    for(index in 1:length(possible)){
      transi <- length(intersect(data[[possible[index]]][[5]], activ[[5]]))
      if(transi>0){
        if(length(data[[possible[index]]])==5){
          data[[possible[index]]][[6]] <- rbind(c(check,transi), NULL)
        } else{
          data[[possible[index]]][[6]] <- rbind(data[[possible[index]]][[6]], c(check,transi) )
        }

      }


    }
  }


  ## Eingaenge

  relevant <- (a[,2]<max(a[,2]))
  relevant <- (1:nrow(a))[relevant]

  for(check in relevant){
    nr <- 1
    poss <- data[[check]][[6]][,1]
    for(index in poss){
      if(length(data[[index]])<7){
        data[[index]][[7]] <- matrix(c(check, data[[check]][[6]][nr,2]), nrow=1)
      } else{
        data[[index]][[7]] <- rbind( data[[index]][[7]], c(check, data[[check]][[6]][nr,2]))
      }

      nr <- nr +1
    }
  }
  return(data)

}

'
