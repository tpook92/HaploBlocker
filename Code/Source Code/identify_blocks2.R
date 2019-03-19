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
#' @export

identify_blocks2 <- function(data, indi, nwindow, mindestanteil=0.95, consider_knoten=TRUE, consider_trans=TRUE,
                            trans_min=5, subgroups=NULL, min_per_subgroup= 0, intersect_func=intersect){
  if(length(subgroups)==0){
    subgroups <- list(1:indi)
  }
  b <- knoten_size(data) # Individuen pro Block
  blocklist <- list()
  blockcounter <- 1
  for(start in 1:length(data)){
    options <- list()
    if(consider_knoten==TRUE){
      for(group in subgroups){
        options[[length(options)+1]] <- list(start, data[[start]][[4]], intersect_func(data[[start]][[5]], group))
      }

    }
    if(consider_trans==TRUE){
      valid <- data[[start]][[6]][(data[[start]][[6]][,2]>=trans_min) * (data[[start]][[6]][,1]>0) * 1:nrow(data[[start]][[6]]),1]
      if(length(valid)>0){
        for(index in valid){
          for(group in subgroups){
            options[[length(options)+1]] <- list(c(start), c(data[[start]][[4]]), intersect_func(group, intersect_func(data[[start]][[5]], data[[index]][[5]]) ))
          }
        }
      }
    }

    if(min_per_subgroup>0){
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
    for(option in options){
      block <- option[[1]]
      windows <- option[[2]]
      individuen <- option[[3]]
      checker <- 1
      last <- start
      # Verlaengerung Nach Hinten
      while(checker==1){
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
          if(last.last==last) checker <- 0
        } else{
          checker <- 0
        }

      }
      # Verlaengerung Nach Vorne
      # individuen <- data[[start]][[5]]
      checker <- 1
      last <- start
      while(checker==1){
        last.last <- last
        if(data[[last]][[1]]$window>1){
          possible <- unique(c(0, data[[last]][[7]][,1]))[-1]
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
              last <- xyz
              individuen <- intersect_func(individuen, data[[xyz]][[5]])
            }
          }
          if(last.last==last) checker <- 0
        } else{
          checker <- 0
        }

      }

      blocklist[[blockcounter]] <- list()
      blocklist[[blockcounter]][[1]] <- block
      blocklist[[blockcounter]][[2]] <- data[[min(block)]][[1]]
      blocklist[[blockcounter]][[3]] <- data[[max(block)]][[2]]
      blocklist[[blockcounter]][[4]] <- windows
      blocklist[[blockcounter]][[6]] <- individuen
      blocklist[[blockcounter]][[5]] <- length(individuen)
      blockcounter <- blockcounter + 1
    }
  }

  return(blocklist)
}
