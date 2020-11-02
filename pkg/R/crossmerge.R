#' Cross Merge nodes
#'
#' Internal Function perform a cross merge for node data
#' @param data1 node-dataset
#' @param indi number of haplotypes in the dataset
#' @param nwindow number of windows in the dataset
#' @param intersect_func Used intersect-function (internally relevant for computation time)
#' @param a Internal helper for a efficient abort criterion
#' @return window cluster after simple-merge

crossmerge <- function(data1, indi, nwindow, a=NULL, intersect_func=intersect){

  if(length(a)==0){
    a <- start_end_block(data1)
  }

  relevant <- (a[,1]>1) * (a[,2]<max(a[,2])) * (1:nrow(a))
  relevant <- (1:nrow(a))[relevant]

  while(length(relevant)>0){
    index <- relevant[1]
    if(data1[[index]][[1]]$window >1 && data1[[index]][[2]]$window <nwindow && data1[[index]][[3]]>0){
      eingang <- data1[[index]][[7]][,1]
      ausgang <- data1[[index]][[6]][,1]
    } else{
      eingang <- ausgang <- NULL
    }
    if(length(eingang)>0 && length(ausgang)>0){
      emerged <- numeric(length(eingang))
      amerged <- numeric(length(ausgang))
      d6 <- data1[[index]][[6]][,2]
      d7 <- data1[[index]][[7]][,2]

      no_change <- TRUE
      for(index2 in 1:length(eingang)){
        for(index3 in 1:length(ausgang)){

          # Check Ob Ausgang/Eingang-Kombination Merge-Kandidat ist
          if(amerged[index3]==0 && emerged[index2]==0){
              activ_a <- ausgang[index3]
              activ_e <- eingang[index2]

              if(no_change){
                pos_a <- index3
                pos_e <- index2
              } else{
                pos_a <- which(data1[[index]][[6]][,1]==activ_a)
                pos_e <- which(data1[[index]][[7]][,1]==activ_e)
              }



              #if(data1[[index]][[6]][pos_a,2]== data1[[index]][[7]][pos_e,2]){
              if(d6[pos_a]== d7[pos_e]){

                # Bestimmung der Exakten Aus-und Eingangstiere
                if(data1[[index]][[7]][pos_e,1]==0){
                  sonstige <- NULL
                  if(pos_e>1){
                    for(temp in 1:(pos_e-1)){
                      sonstige <- c(sonstige, data1[[data1[[index]][[7]][temp,1]]][[5]])
                    }
                  }
                  if(length(sonstige)>0){
                    etier <- intersect_func((1:indi)[-sonstige], data1[[index]][[5]])
                  } else{
                    etier <- data1[[index]][[5]]
                  }

                } else{
                  etier <- intersect_func(data1[[data1[[index]][[7]][pos_e,1]]][[5]], data1[[index]][[5]])
                }
                if(data1[[index]][[6]][pos_a,1]==0){
                  sonstige <- NULL
                  if(pos_a>1){
                    for(temp in 1:(pos_a-1)){
                      sonstige <- c(sonstige, data1[[data1[[index]][[6]][temp,1]]][[5]])
                    }
                  }
                  if(length(sonstige)>0){
                    atier <- intersect_func((1:indi)[-sonstige], data1[[index]][[5]])
                  } else{
                    atier <- data1[[index]][[5]]
                  }

                } else{
                  atier <- intersect_func(data1[[data1[[index]][[6]][pos_a,1]]][[5]], data1[[index]][[5]])
                }

                # Merge Knoten Wenn Aus/Eingangstiere Identisch
                if(length(intersect_func(etier, atier))==max(length(atier), length(etier))){
                  no_change <- FALSE
                  if(activ_a==0 && activ_e==0){
                    # 1. Fall: Eingangsknoten 0, Ausgangsknoten 0 - Erstelle zusaetzlichen Knoten
                    new_k <- length(data1) +1
                    data1[[new_k]] <- data1[[index]]
                    data1[[new_k]][[5]] <- atier
                    data1[[index]][[5]] <- intersect_func(data1[[index]][[5]], (1:indi)[-atier])
                    data1[[new_k]][[3]] <- length(data1[[new_k]][[5]])
                    data1[[index]][[3]] <- length(data1[[index]][[5]])
                    data1[[new_k]][[6]] <- cbind(0,data1[[new_k]][[3]])
                    data1[[new_k]][[7]] <- cbind(0,data1[[new_k]][[3]])

                    data1[[index]] <- calculate_new_transition(data1, index, nwindow, intersect_func=intersect_func)
                    amerged[index3] <- 1
                    emerged[index2] <- 1
                  } else if(activ_a==0){
                    # 2. Fall Eingangsknoten Existiert, Ausgangsknoten 0
                    # Merge nur wenn auch alle Tiere im Eingangsknoten nach Index uebergehen.
                    if(data1[[activ_e]][[3]]==length(etier)){
                      data1[[activ_e]][[2]] <- data1[[index]][[2]]
                      data1[[activ_e]][[4]] <- c(data1[[activ_e]][[4]], data1[[index]][[4]])
                      data1[[activ_e]][[6]] <- cbind(0,data1[[activ_e]][[3]])

                      data1[[index]][[5]] <- intersect_func(data1[[index]][[5]], (1:indi)[-etier])
                      data1[[index]][[3]] <- length(data1[[index]][[5]])
                      data1[[index]] <- calculate_new_transition(data1, index, nwindow, intersect_func=intersect_func)
                      amerged[index3] <- 1
                      emerged[index2] <- 1
                    }
                  } else if(activ_e==0){
                    # 3. Fall Eingangsknoten existiert Nicht, Ausgangsknoten existiert und enthaelt Exakt knoten[index]
                    # Merge analog zur 2.Fall nur wen alle Tiere uebergehen.
                    if(length(atier)== data1[[activ_a]][[3]]){
                      data1[[activ_a]][[1]] <- data1[[index]][[1]]
                      data1[[activ_a]][[4]] <- c(data1[[index]][[4]], data1[[activ_a]][[4]])
                      data1[[activ_a]][[7]] <- cbind(0,data1[[activ_a]][[3]])

                      data1[[index]][[5]] <- intersect_func(data1[[index]][[5]], (1:indi)[-etier])
                      data1[[index]][[3]] <- length(data1[[index]][[5]])
                      data1[[index]] <- calculate_new_transition(data1, index, nwindow, intersect_func=intersect_func)
                      amerged[index3] <- 1
                      emerged[index2] <- 1
                    }
                  } else if(length(atier)== data1[[activ_a]][[3]]){
                    # 4. Fall Ein- und Ausgangsknoten existieren
                    # saemtliche Tiere Aus Ausgangsknoten kommen aus knoten[index]
                    data1[[index]][[5]] <- intersect_func(data1[[index]][[5]], (1:indi)[-etier])
                    data1[[activ_e]][[5]] <- intersect_func(data1[[activ_e]][[5]], (1:indi)[-etier])
                    data1[[index]][[3]] <- length(data1[[index]][[5]])
                    data1[[activ_e]][[3]] <- length(data1[[activ_e]][[5]])
                    data1[[activ_a]][[1]] <- data1[[activ_e]][[1]]
                    data1[[activ_a]][[4]] <- c(data1[[activ_e]][[4]], data1[[index]][[4]], data1[[activ_a]][[4]])
                    if(data1[[activ_a]][[1]]$window>1){
                      data1[[activ_a]][[7]] <- data1[[activ_e]][[7]] # Enthaelt alle moeglichen Eingaenge von activ_a und mehr!
                    } else{
                      data1[[activ_a]][[7]] <- cbind(0,0)
                    }
                    data1[[index]] <- calculate_new_transition(data1,index,nwindow, intersect_func=intersect_func)
                    data1[[activ_e]] <- calculate_new_transition(data1, activ_e, nwindow, intersect_func=intersect_func)
                    data1[[activ_a]] <- calculate_new_transition(data1, activ_a, nwindow, intersect_func=intersect_func)
                    amerged[index3] <- 1
                    emerged[index2] <- 1
                    # samtliche Knoten die knoten[activ_e] gehen nun in knoten[activ_a]
                    recalc <- data1[[activ_a]][[7]][,1]
                    if(length(recalc[recalc>0])>0){
                      for(rindex in recalc[recalc>0]){
                        data1[[rindex]] <- calculate_new_transition(data1, rindex, nwindow, add_a = activ_a, intersect_func=intersect_func)
                      }
                    }
                  } else if(length(atier)== data1[[activ_e]][[3]]){
                    # 5. Fall Ein- und Ausgangsknoten existieren
                    # saemtliche Tiere im Eingangsknoten gehen in knoten[index]
                    data1[[index]][[5]] <- intersect_func(data1[[index]][[5]], (1:indi)[-etier])
                    data1[[activ_a]][[5]] <- intersect_func(data1[[activ_a]][[5]], (1:indi)[-etier])
                    data1[[index]][[3]] <- length(data1[[index]][[5]])
                    data1[[activ_a]][[3]] <- length(data1[[activ_a]][[5]])
                    data1[[activ_e]][[2]] <- data1[[activ_a]][[2]]
                    data1[[activ_e]][[4]] <- c(data1[[activ_e]][[4]], data1[[index]][[4]], data1[[activ_a]][[4]])
                    data1[[activ_e]][[6]] <- data1[[activ_a]][[6]] # Enthaelt alle moeglichen Eingaenge von activ_a und mehr!
                    data1[[index]] <- calculate_new_transition(data1,index,nwindow, intersect_func=intersect_func)
                    data1[[activ_e]] <- calculate_new_transition(data1, activ_e, nwindow, intersect_func=intersect_func)
                    data1[[activ_a]] <- calculate_new_transition(data1, activ_a, nwindow, intersect_func=intersect_func)
                    amerged[index3] <- 1
                    emerged[index2] <- 1

                    # samtliche Knoten die knoten[activ_a] gehen nun in knoten[activ_e]
                    recalc <- data1[[activ_e]][[6]][,1]
                    if(length(recalc[recalc>0])>0){
                      for(rindex in recalc[recalc>0]){
                        data1[[rindex]] <- calculate_new_transition(data1, rindex, nwindow, add_e = activ_e, intersect_func=intersect_func)
                      }
                    }

                  }





                }
              }
          }
        }
      }


    }
    relevant <- relevant[-1]
  }
  data1 <- renaming_combi(data1, nwindow)
}
