#' Simple Merge nodes
#'
#' Internal Function perform a simple merge for node data
#' @param data node-dataset
#' @param indi number of haplotypes in the dataset
#' @param nwindow number of windows in the dataset
#' @param intersect_func Used intersect-function (internally relevant for computation time)
#' @export


simple_merge_prob <- function(data, indi, nwindow, intersect_func=intersect){

  a <- start_end_block(data)
  relevant <- (a[,1]>1)
  relevant <- (1:nrow(a))[relevant]
  relevant <- relevant[length(relevant):1]
  for(check in relevant){
    activ <- data[[check]]
    #eingaenge <- unique(c(0,activ[[7]][,1]))[-1]
    eingaenge <- activ[[7]][,1]
    if(length(eingaenge)==1 && eingaenge!=0){
      #ausgaenge <- unique(c(0,data[[eingaenge]][[6]][,1]))[-1]
      ausgaenge <- data[[eingaenge]][[6]][,1]
      if(length(ausgaenge)==1 && ausgaenge!=0){
        if(data[[eingaenge]][[3]]==activ[[3]]){
          old <- data[[eingaenge]]
          # Merge-Prozess!
          data[[eingaenge]][[2]] <- activ[[2]]
          data[[eingaenge]][[4]] <- c(old[[4]], activ[[4]])
          data[[eingaenge]][[5]] <- intersect_func(old[[5]], activ[[5]])
          data[[eingaenge]][[3]] <- length(data[[eingaenge]][[5]])

          # Wahrscheinlichkeiten Anpassung

          ausgaenge_old <- unique(c(0,activ[[6]][,1]))[-1]
          data[[check]][[5]] <- numeric(0)
          data[[check]][[3]] <- 0

          if(data[[eingaenge]][[2]]$window<nwindow){
            data[[eingaenge]] <- calculate_new_transition(data, eingaenge, nwindow, add_a = activ[[6]][,1], intersect_func=intersect_func)
          } else{
            data[[eingaenge]] <- calculate_new_transition(data, eingaenge, nwindow, intersect_func=intersect_func)
          }

          for(switch in ausgaenge_old){
            data[[switch]] <- calculate_new_transition(data, switch, nwindow, add_e = eingaenge, intersect_func=intersect_func )
          }

        }
      }
    }
  }
  data <- renaming_combi(data, nwindow)
  return(data)

}
