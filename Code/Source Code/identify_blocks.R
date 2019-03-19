#' Identify possible blocks base on the node-dataset
#'
#' Function to identify possible blocks base on the node-dataset
#' @param data node-dataset
#' @param indi number of haplotypes in the dataset
#' @param nwindow number of windows in the dataset
#' @param min_share minimum percentage of transition to the same block for extension (default: 0.95, step III)
#' @param consider_nodes Use nodes to identify blocks (default: TRUE)
#' @param consider_edge Use edges between nodes to identify blocks (default: TRUE)
#' @param edge_min minimum number of haplotypes per transition to use in consider_edge (default: 5)
#' @param subgroups possible subgroups to consider in the block identification process (default: NULL - list(1:indi))
#' @param min_per_subgroup minimum number of haplotypes per block per subgroup (default: 0)
#' @param subgroup_exception allow for not all subgroups to include a block
#' @param intersect_func Used intersect-function (internally relevant for computation time)
#' @param identify_filter If TRUE filter out Starting blocks to avoid duplicated (default: FALSE)
#' @param consider_multi If TRUE use multi-level edges to identify blocks (default: FALSE)
#' @param multi_min minimum number of haplotypes per multi transition to use in consider_multi (default: 5)
#' @param double_share Extended-block-identification share t
#' @param node_min minimum number of haplotypes per block (default: 5)
#' @export

identify_blocks <- function(data, indi, nwindow, min_share=0.95, consider_nodes=TRUE, consider_edge=TRUE,
                            edge_min=5, subgroups=NULL, min_per_subgroup= 0, intersect_func=intersect,
                            identify_filter=FALSE, consider_multi=FALSE,
                            multi_min=5, subgroup_exception=0,
                            double_share=1, node_min=5){
  if(length(subgroups)==0){
    subgroups <- list(1:indi)
  }
  b <- nodes_size(data) # Individuen pro Block
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
    if(consider_nodes==TRUE){
      if(nrow(data[[start]][[7]])>1 || (nrow(data[[start]][[7]])==1 && data[[start]][[7]][1,1]==0) || (nrow(data[[start]][[7]])==1 && nrow(data[[data[[start]][[7]][1,1]]][[6]])>=2)){
        for(group in identify_subgroups){
          options[[length(options)+1]] <- list(start, data[[start]][[4]], intersect_func(data[[start]][[5]], group), start, start)
        }
      }
    }
    if(consider_edge==TRUE){
      valid <- data[[start]][[6]][(data[[start]][[6]][,2]>=edge_min) * (data[[start]][[6]][,1]>0) * 1:nrow(data[[start]][[6]]),1]
      if((nrow(data[[start]][[6]])!=1 && data[[start]][[6]][1,1]!=0 ) || consider_nodes==FALSE){
        if(length(valid)>0){
          for(index in valid){
            for(group in identify_subgroups){
              options[[length(options)+1]] <- list(c(start), c(data[[start]][[4]]), intersect_func(group, intersect_func(data[[start]][[5]], data[[index]][[5]])), start, start )
            }
          }
        }
      }
    }
    if(consider_multi==TRUE){
      valid <- data[[start]][[6]][(data[[start]][[6]][,2]>=edge_min) * (data[[start]][[6]][,1]>0) * 1:nrow(data[[start]][[6]]),1]
      if((nrow(data[[start]][[6]])!=1 && data[[start]][[6]][1,1]!=0 ) || consider_nodes==FALSE){
        if(length(valid)>0){
          for(index in valid){
            valid2 <- data[[index]][[6]][(data[[index]][[6]][,2]>=edge_min) * (data[[index]][[6]][,1]>0) * 1:nrow(data[[index]][[6]]),1]

            if((nrow(data[[index]][[6]])!=1 && data[[index]][[6]][1,1]!=0 ) || consider_nodes==FALSE){
              if(length(valid2)>0){
                for(index2 in valid2){
                  for(group in identify_subgroups){
                    options[[length(options)+1]] <- list(c(start), c(data[[start]][[4]]), intersect_func(intersect_func(group, intersect_func(data[[start]][[5]], data[[index]][[5]])), data[[index2]][[5]]), start, start)
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
          exception <- subgroup_exception
          for(group in subgroups){
            if(min_per_subgroup> length(intersect_func(group, options[[index]][[3]]))){
              exception <- exception -1
              if(exception<0){
                checker <- 0
              }
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

    new_options <- options
    if(start>1){
      se <- blocklist_startend(blocklist, type = "window")
      size <- blocklist_size(blocklist)
    }
    new_blocks <- list()
    new_blocks_counter <- 1
    while(length(new_options)>0){
      options <- unique(new_options)
      new_options <- list()

      # option validation
      if(length(options[[1]])>=5 && start>1){

        for(index in length(options):1){
          size_check <- size == length(options[[index]][[3]])
          se_check1 <- se[,1] <= data[[options[[index]][[4]]]][[2]][[1]]
          se_check2 <- se[,2] >= data[[options[[index]][[5]]]][[2]][[1]]
          candidates <- which((size_check * se_check1 * se_check2)==1)
          ncan <- length(candidates)
          index2 <- 1
          no.merge <- 1
          while(no.merge==1 && index2 <= ncan){
            if(length(intersect_func(options[[index]][[3]], blocklist[[candidates[index2]]][[6]])) == length(options[[index]][[3]])){
              no.merge <- 0
              options[[index]] <- NULL
            } else{
              index2 <- index2 + 1
            }
          }
        }

      }

      for(option in options){
        block <- option[[1]]
        windows <- option[[2]]
        individuen <- option[[3]]
        checker1 <- 1
        checker2 <- 1

        if(length(option)>=5){
          last <- option[[5]]
          first <- option[[4]]
        } else{
          last <- start
          first <- start
        }


        while((checker1+checker2)>0){
          if(checker1==1){
            last.last <- last
            if(data[[last]][[2]]$window <nwindow){
              possible <- unique(c(0, data[[last]][[6]][,1]))[-1]
              for(xyz in possible){
                anteil <- length(intersect_func(individuen, data[[xyz]][[5]]))/length(individuen)
                subgroup_stop <- 1
                if(min_per_subgroup>0){
                  exception <- subgroup_exception
                  for(groups in subgroups){
                    anzahl <- length(intersect_func(intersect_func(individuen,data[[xyz]][[5]] ), groups))
                    if(anzahl<min_per_subgroup){
                      exception <- exception -1
                      if(exception<0){
                        subgroup_stop <- 0
                      }
                    }
                  }
                }
                if(subgroup_stop==1 && anteil >= min_share){
                  block <- c(block, xyz)
                  windows <- c(windows, data[[xyz]][[4]])
                  last <- xyz
                  individuen <- intersect_func(individuen, data[[xyz]][[5]])
                } else if(subgroup_stop==1 && anteil >= double_share && length(intersect_func(individuen, data[[xyz]][[5]]))>=node_min){
                  tempn <- length(new_options)+1
                  new_options[[tempn]] <- list()
                  new_options[[tempn]][[1]] <- c(block, xyz)
                  new_options[[tempn]][[2]] <- c(windows, data[[xyz]][[4]])
                  new_options[[tempn]][[3]] <- intersect_func(individuen, data[[xyz]][[5]])
                  new_options[[tempn]][[4]] <- first
                  new_options[[tempn]][[5]] <- xyz
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
                  exception <- subgroup_exception
                  for(groups in subgroups){
                    anzahl <- length(intersect_func(intersect_func(individuen,data[[xyz]][[5]] ), groups))
                    if(anzahl<min_per_subgroup){
                      exception <- exception -1
                      if(exception<0){
                        subgroup_stop <- 0
                      }
                    }
                  }
                }
                if(subgroup_stop==1 && anteil >= min_share){
                  block <- c(xyz,block)
                  windows <- c(data[[xyz]][[4]], windows)
                  first <- xyz
                  individuen <- intersect_func(individuen, data[[xyz]][[5]])
                } else if(subgroup_stop==1 && anteil >= double_share && length(intersect_func(individuen, data[[xyz]][[5]]))>=node_min){
                  tempn <- length(new_options)+1
                  new_options[[tempn]] <- list()
                  new_options[[tempn]][[1]] <- c(xyz,block)
                  new_options[[tempn]][[2]] <- c(data[[xyz]][[4]], windows)
                  new_options[[tempn]][[3]] <- intersect_func(individuen, data[[xyz]][[5]])
                  new_options[[tempn]][[4]] <- xyz
                  new_options[[tempn]][[5]] <- last
                }
              }
              if(last.first==first) checker2 <- 0
            } else{
              checker2 <- 0
            }

          }

        }
        # Verlaengerung Nach Hinten

        new_blocks[[new_blocks_counter]] <- list()
        new_blocks[[new_blocks_counter]][[1]] <- block
        new_blocks[[new_blocks_counter]][[2]] <- data[[min(block)]][[1]]
        new_blocks[[new_blocks_counter]][[3]] <- data[[max(block)]][[2]]
        new_blocks[[new_blocks_counter]][[4]] <- windows
        new_blocks[[new_blocks_counter]][[6]] <- individuen
        new_blocks[[new_blocks_counter]][[5]] <- length(individuen)
        new_blocks_counter <- new_blocks_counter + 1


      }



    }

    new_blocks <- unique(new_blocks)

    if(length(new_blocks)>0){
      for(index3 in 1:length(new_blocks)){
        blocklist[[blockcounter]] <- list()
        blocklist[[blockcounter]][[1]] <- new_blocks[[index3]][[1]]
        blocklist[[blockcounter]][[2]] <- new_blocks[[index3]][[2]]
        blocklist[[blockcounter]][[3]] <- new_blocks[[index3]][[3]]
        blocklist[[blockcounter]][[4]] <- new_blocks[[index3]][[4]]
        blocklist[[blockcounter]][[6]] <- new_blocks[[index3]][[6]]
        blocklist[[blockcounter]][[5]] <- new_blocks[[index3]][[5]]
        blockcounter <- blockcounter + 1
      }
    }

    ## FIND THE ERROR!
    #blocklist <- unique(blocklist)
    #blockcounter <- length(blocklist)+1

    if(identify_filter){
      activ_blocks <- activ_blocks[activ_blocks[,2]>=data[[start]]$end$window,]
      if(length(options)>0){
        activ_blocks <- rbind(activ_blocks, cbind(blocklist_startend(blocklist, first_block = length(blocklist)- length(options)+1), (length(blocklist)- length(options)+1):length(blocklist), blocklist_size(blocklist, first_block = length(blocklist)- length(options)+1)))
      }
    }
  }

  return(blocklist)
}
