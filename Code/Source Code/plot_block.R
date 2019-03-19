#' Plot Block_structure
#'
#' Function to create a graphical representation of the blocklist
#' @param oriantation Position on with to sort haplotypes for similarity (default: "snp", options: "front", "mid","back")
#' @param snp_ori Select SNP for the oriantation algorithm if oriantation=="snp"
#' @param export_order Export the chosen haplotype order after sorting
#' @param import_order Import a given haplotype order
#' @param min_to_plot Only include blocks with at least this many haplotypes in it (default:5)
#' @param intensity Color-intensity for the plotted blocks (default: 0.5 (options: 0-1))
#' @param include Include haplotypes in no block in the oriantation position in the graph
#' @param indi number of haplotypes in the dataset
#' @param type length measure (default: "window" , alt: "snp", "bp")
#' @param max_step Maximum distance of blocks to consider in the oriantation step (default: 500)
#' @param blocklist block-dataset
#' @param add_sort If FALSE deactivate sorting haplotype for local similarity
#' @export


plot_block <- function(blocklist, type="snp", oriantation="snp", include=TRUE, indi=313, min_to_plot = 5,
                      intensity=0.5, add_sort=TRUE, max_step=500,
                      snp_ori=NULL,
                      export_order=FALSE, import_order=FALSE){
  oriantation_f <- oriantation
  se <- blocklist_startend(blocklist, type=type)
  order <- NULL

  if(length(oriantation_f)==0){
    order <- 1:indi
  } else if(oriantation_f[1]=="front"){
    oriantation <- 1:8
  } else if(oriantation_f[1]=="mid"){
    oriantation <- 1:8 + ceiling(nrow(se)/2-3)
  } else if(oriantation_f[1]=="back"){
    oriantation <- 1:8 + nrow(se) -8
  } else if(oriantation_f[1]=="snp"){
    if(length(snp_ori)==0){
      snp_ori <- mean(c(max(se), min(se)))
    }
    oriantation <- which(((se[,1]<= snp_ori)+ (se[,2]>= snp_ori))==2)
  } else{
    order <- oriantation_f
  }

  if(length(oriantation)>0 && import_order==FALSE){
    sorted <- oriantation[1]
    oriantation <- oriantation[-1]
    while(length(oriantation)>0){
      overlap <- numeric(length(oriantation))
      actives <- blocklist[[sorted[length(sorted)]]][[6]]
      i1 <- 1
      for(index in oriantation){
        nextone <- blocklist[[index]][[6]]
        overlap[i1] <- length(intersect(nextone,actives))/ length(nextone)
        i1 <- i1 + 1
      }
      sorted <- c(sorted, oriantation[which.max(overlap)[1]])
      oriantation <- oriantation[-which.max(overlap)[1]]
    }
    oriantation <- sorted

  }
  if(length(order)==0){
    if(include==TRUE){
      oriantation <- c(oriantation, 0)
    }
    for(index in oriantation){
      print(index)
      pl <- length(order)
      if(index==0){
        order <- unique(c(order,1:indi))
        index <- ceiling(median(oriantation[-length(oriantation)]))
      } else{
        order <- unique(c(order,blocklist[[index]][[6]]))
      }

      if(add_sort==TRUE && length(order)>pl){
        if(pl>0){
          added <- order[-(1:pl)]
        } else{
          added <- order
        }
        groups <- list()
        groups[[1]] <- added
        for(step in 1:max_step){

          if(length(groups)<length(added)){
            if((index-step) >0 && length(intersect(blocklist[[index-step]][[6]], added))>0){
              new_groups <- list()
              for(group in 1:length(groups)){
                same <- base::intersect(groups[[group]], blocklist[[index-step]][[6]])
                rest <- base::intersect(groups[[group]], (1:indi)[-blocklist[[index-step]][[6]]])
                if(length(same)>0){
                  new_groups[[length(new_groups)+1]] <- same
                }
                if(length(rest)>0){
                  new_groups[[length(new_groups)+1]] <- rest
                }

              }
              groups <- new_groups
            }
            if((index+step) <=length(blocklist) && length(intersect(blocklist[[index+step]][[6]], added))>0){
              new_groups <- list()
              for(group in 1:length(groups)){
                same <- base::intersect(groups[[group]], blocklist[[index+step]][[6]])
                rest <- base::intersect(groups[[group]], (1:indi)[-blocklist[[index+step]][[6]]])
                if(length(same)>0){
                  new_groups[[length(new_groups)+1]] <- same
                }
                if(length(rest)>0){
                  new_groups[[length(new_groups)+1]] <- rest
                }


              }
              groups <- new_groups
            }
          }

        }
        if(pl>0){
          order <- c(order[(1:pl)],unlist(groups))
        } else{
          order <- c(unlist(groups))
        }

        print(unlist(groups))
        print(order)

      }
    }


  }
  plot(0,-1000,ylim=c(0,length(order)), xlim=c(1,max(se)), ylab="haplotype", xlab="SNP",
       cex.axis=1.3, cex.lab=1.3)

  activ_colors <- numeric(length(blocklist))
  for(index in 1:length(blocklist)){
    overlap <- duplicated(c(blocklist[[index]][[6]], order))[-(1:blocklist[[index]][[5]])]
    if(sum(overlap) >= min_to_plot){
      print(index)
      taken <- base::intersect(which(se[,1]<=se[index,2]), which(se[,2]>=se[index,1]))
      block_color <- sort(unique(c(0,activ_colors[taken])))
      activ_colors[index] <- min(length(block_color),which(block_color!=(1:length(block_color)-1))-1)
      for(index2 in which(overlap)){
        polygon(c(se[index,1], se[index,2], se[index,2], se[index,1]), index2-c(1,1,0,0),
                col=adjustcolor(activ_colors[index],alpha.f=intensity), lty=0)
      }
    }
  }
  if(export_order==TRUE){
    return(order)
  }
}
