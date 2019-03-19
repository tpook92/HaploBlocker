#' Simple Merge knots
#'
#' Internal Function to calculate transition probabilities for a single knot
#' @param data knot-dataset
#' @param knoten knot for which to calculate transition probabilities
#' @param nwindow number of windows in the dataset
#' @param add_e additional possible previous knots
#' @param add_a additional possible following knots
#' @param intersect_func Used intersect-function (internally relevant for computation time)
#' @export

calculate_new_transition <- function(data, knoten, nwindow, add_e= NULL, add_a= NULL, intersect_func=intersect){
  new.knoten <- data[[knoten]]
  if(new.knoten[[1]]$window>1 && new.knoten[[3]]>0){
    eingaenge <- sort(unique(c(0, new.knoten[[7]][(new.knoten[[7]][,1]!=0),1], add_e))[-1], decreasing=TRUE)
    haeufig <- numeric(length(eingaenge))
    if(length(haeufig)>0){
      for(index in 1:length(haeufig)){
        haeufig[index] <- length(intersect_func(data[[eingaenge[index]]][[5]], new.knoten[[5]]))
      }
    }
    eingaenge <- eingaenge[haeufig>0]
    haeufig <- haeufig[haeufig>0]
    if(sum(haeufig)< new.knoten[[3]]){
      new.knoten[[7]] <- rbind(cbind(eingaenge, haeufig), c(0, new.knoten[[3]]- sum(haeufig)))
    } else{
      new.knoten[[7]] <- cbind(eingaenge, haeufig)
    }
  } else{
    new.knoten[[7]] <- cbind(0,0)
  }
  if(new.knoten[[2]]$window <nwindow && new.knoten[[3]]>0){
    ausgaenge <- sort(unique(c(0, new.knoten[[6]][(new.knoten[[6]][,1]!=0),1], add_a))[-1], decreasing=TRUE)
    haeufig <- numeric(length(ausgaenge))
    if(length(haeufig)>0){
      for(index in 1:length(haeufig)){
        haeufig[index] <- length(intersect_func(data[[ausgaenge[index]]][[5]], new.knoten[[5]]))
      }
    }
    ausgaenge <- ausgaenge[haeufig>0]
    haeufig <- haeufig[haeufig>0]
    if(sum(haeufig)< new.knoten[[3]]){
      new.knoten[[6]] <- rbind(cbind(ausgaenge, haeufig), c(0, new.knoten[[3]]- sum(haeufig)))
    } else{
      new.knoten[[6]] <- cbind(ausgaenge, haeufig)
    }
  } else{
    new.knoten[[6]] <- cbind(0,0)
  }
  return(new.knoten)
}
