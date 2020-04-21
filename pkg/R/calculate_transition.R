#' Calculate transition probabilities between nodes
#'
#' Internal Function to calculate transition probabilities between nodes
#' @param data node-dataset
#' @param intersect_func Used intersect-function (internally relevant for computation time)
#' @return window cluster


calculate_transition <- function(data, intersect_func=intersect){

  a <- start_end_block(data)
  relevant <- (a[,1]>1)
  relevant <- (1:nrow(a))[relevant]
  relevant <- relevant[length(relevant):1]
  for(index in which(a[,2]==max(a[,2]))){
    data[[index]][[6]] <- cbind(0,0)
  }
  for(check in relevant){
    activ <- data[[check]]
    possible <- (1:nrow(a))[(a[,2]==activ[[1]]$window-1)]

    for(index in 1:length(possible)){
      transi <- length(intersect(data[[possible[index]]][[5]], activ[[5]]))
      if(transi>0){
        if(length(data[[possible[index]]])==5){
          data[[possible[index]]][[6]] <- cbind(check,transi)
        } else{
          data[[possible[index]]][[6]] <- rbind(data[[possible[index]]][[6]], c(check,transi) )
        }

      }


    }
  }


  ## Eingaenge

  relevant <- (a[,2]<max(a[,2]))
  relevant <- (1:nrow(a))[relevant]

  for(index in which(a[,1]==1)){
    data[[index]][[7]] <- cbind(0,0)
  }

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
