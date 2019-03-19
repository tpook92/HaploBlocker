#' Identify possible blocks base on the knot-dataset
#'
#' Function to identify possible blocks base on the knot-dataset
#' @param data knot-dataset
#' @param indi number of haplotypes in the dataset
#' @param nwindow number of windows in the dataset
#' @param mindestanteil minimum percentage of transition to the same block for extension (default: 0.95, step III)
#' @param consider_knoten Use knots to identify blocks (default: TRUE)
#' @param consider_trans Use edges between knots to identify blocks (default: TRUE)
#' @param trans_min minimum number of haplotypes per transition to use in consider_trans (default: 5)
#' @param subgroups possible subgroups to consider in the block identification process (default: NULL - list(1:indi))
#' @param min_per_subgroup minimum number of haplotypes per block per subgroup (default: 0)
#' @param intersect_func Used intersect-function (internally relevant for computation time)
#' @param identify_filter If TRUE filter out Starting blocks to avoid duplicated (default: FALSE)
#' @param consider_multi If TRUE use multi-level edges to identify blocks (default: FALSE)
#' @param multi_min minimum number of haplotypes per multi transition to use in consider_multi (default: 5)
#' @export

identify_blocks3 <- function(data, indi, nwindow, mindestanteil=0.95, consider_knoten=TRUE, consider_trans=TRUE,
                            trans_min=5, subgroups=NULL, min_per_subgroup= 0, intersect_func=intersect,
                            identify_filter=FALSE, consider_multi=FALSE,
                            multi_min=5){
  if(length(subgroups)==0){
    subgroups <- list(1:indi)
  }
  b <- knoten_size(data) # Individuen pro Block
  blocklist <- list()
  blockcounter <- 1
  activ_blocks <- NULL
  if(min_per_subgroup>0){
    identify_subgroups <- list(sort(unique(unlist(subgroups))))
  } else{
    identify_subgroups <- subgroups
  }
  for(start in 1:length(data)){
    options <- list()
    if(consider_knoten==TRUE){
      if(nrow(data[[start]][[7]])>1 || (nrow(data[[start]][[7]])==1 && data[[start]][[7]][1,1]==0) || (nrow(data[[start]][[7]])==1 && nrow(data[[data[[start]][[7]][1,1]]][[6]])>=2)){
        for(group in identify_subgroups){
          options[[length(options)+1]] <- list(start, data[[start]][[4]], intersect_func(data[[start]][[5]], group))
        }
      }
    }
    if(consider_trans==TRUE){
      valid <- data[[start]][[6]][(data[[start]][[6]][,2]>=trans_min) * (data[[start]][[6]][,1]>0) * 1:nrow(data[[start]][[6]]),1]
      if((nrow(data[[start]][[6]])!=1 && data[[start]][[6]][1,1]!=0 ) || consider_knoten==FALSE){
        if(length(valid)>0){
          for(index in valid){
            for(group in identify_subgroups){
              options[[length(options)+1]] <- list(c(start), c(data[[start]][[4]]), intersect_func(group, intersect_func(data[[start]][[5]], data[[index]][[5]])) )
            }
          }
        }
      }
    }
    if(consider_multi==TRUE){
      valid <- data[[start]][[6]][(data[[start]][[6]][,2]>=trans_min) * (data[[start]][[6]][,1]>0) * 1:nrow(data[[start]][[6]]),1]
      if((nrow(data[[start]][[6]])!=1 && data[[start]][[6]][1,1]!=0 ) || consider_knoten==FALSE){
        if(length(valid)>0){
          for(index in valid){
            valid2 <- data[[index]][[6]][(data[[index]][[6]][,2]>=trans_min) * (data[[index]][[6]][,1]>0) * 1:nrow(data[[index]][[6]]),1]

            if((nrow(data[[index]][[6]])!=1 && data[[index]][[6]][1,1]!=0 ) || consider_knoten==FALSE){
              if(length(valid2)>0){
                for(index2 in valid2){
                  for(group in identify_subgroups){
                    options[[length(options)+1]] <- list(c(start), c(data[[start]][[4]]), intersect_func(intersect_func(group, intersect_func(data[[start]][[5]], data[[index]][[5]])), data[[index2]][[5]]))
                    if(length(options[[length(options)]][[3]])<multi_min){
                      options[[length(options)]] <- NULL
                    }
                  }
                }
              }
            }
          }
        }
      }

    }

    if(min_per_subgroup>0){
      if(length(options)>0){
        for(index in length(options):1){
          checker <- 1
          for(group in subgroups){
            if(min_per_subgroup> length(intersect_func(group, options[[index]][[3]]))){
              checker <- 0
            }
          }
          if(checker==0){
            options[[index]] <- NULL
          }
        }
      }
    } else{
      for(index in length(options):1){
        if(length(options[[index]][[3]])==0){
          options[[index]] <- NULL
        }
      }
    }

    if(identify_filter){
      if(length(options)>0){
        for(index in length(options):1){
          candidates <- which(activ_blocks[,4]==length(options[[index]][[3]]))
          for(check in candidates){
            if(prod(options[[index]][[3]]==blocklist[[activ_blocks[check,3]]][[6]])){
              options[[index]] <- NULL
              break

            }
          }
        }
      }

      # Hier werden optionen aussortiert!
    }

    for(option in options){
      block <- option[[1]]
      windows <- option[[2]]
      individuen <- option[[3]]
      checker1 <- 1
      checker2 <- 1

      last <- start
      first <- start

      while((checker1+checker2)>0){
        if(checker1==1){
          last.last <- last
          if(data[[last]][[2]]$window <nwindow){
            possible <- unique(c(0, data[[last]][[6]][,1]))[-1]
            for(xyz in possible){
              anteil <- length(intersect_func(individuen, data[[xyz]][[5]]))/length(individuen)
              subgroup_stop <- 1
              if(min_per_subgroup>0){
                for(groups in subgroups){
                  anzahl <- length(intersect_func(intersect_func(individuen,data[[xyz]][[5]] ), groups))
                  if(anzahl<min_per_subgroup){
                    subgroup_stop <- 0
                  }
                }
              }
              if(subgroup_stop==1 && anteil >= mindestanteil){
                block <- c(block, xyz)
                windows <- c(windows, data[[xyz]][[4]])
                last <- xyz
                individuen <- intersect_func(individuen, data[[xyz]][[5]])
              }
            }
            if(last.last==last) checker1 <- 0
          } else{
            checker1 <- 0
          }

        }

        if(checker2==1){
          last.first <- first
          if(data[[first]][[1]]$window>1){
            possible <- unique(c(0, data[[first]][[7]][,1]))[-1]
            for(xyz in possible){
              anteil <- length(intersect_func(individuen, data[[xyz]][[5]]))/length(individuen)
              subgroup_stop <- 1
              if(min_per_subgroup>0){
                for(groups in subgroups){
                  anzahl <- length(intersect_func(intersect_func(individuen,data[[xyz]][[5]] ), groups))
                  if(anzahl<min_per_subgroup){
                    subgroup_stop <- 0
                  }
                }
              }
              if(subgroup_stop==1 && anteil >= mindestanteil){
                block <- c(xyz,block)
                windows <- c(data[[xyz]][[4]], windows)
                first <- xyz
                individuen <- intersect_func(individuen, data[[xyz]][[5]])
              }
            }
            if(last.first==first) checker2 <- 0
          } else{
            checker2 <- 0
          }

        }

      }
      # Verlaengerung Nach Hinten



      blocklist[[blockcounter]] <- list()
      blocklist[[blockcounter]][[1]] <- block
      blocklist[[blockcounter]][[2]] <- data[[min(block)]][[1]]
      blocklist[[blockcounter]][[3]] <- data[[max(block)]][[2]]
      blocklist[[blockcounter]][[4]] <- windows
      blocklist[[blockcounter]][[6]] <- individuen
      blocklist[[blockcounter]][[5]] <- length(individuen)
      blockcounter <- blockcounter + 1
    }

    if(identify_filter){
      activ_blocks <- activ_blocks[activ_blocks[,2]>=data[[start]]$end$window,]
      if(length(options)>0){
        activ_blocks <- rbind(activ_blocks, cbind(blocklist_startend(blocklist, first_block = length(blocklist)- length(options)+1), (length(blocklist)- length(options)+1):length(blocklist), blocklist_size(blocklist, first_block = length(blocklist)- length(options)+1)))
      }
    }
  }

  return(blocklist)
}
