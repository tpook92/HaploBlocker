#' Plot Block_structure
#'
#' Function to create a graphical representation of the blocklist
#' @param orientation Position on with to sort haplotypes for similarity (default: "snp", options: "front", "mid","back")
#' @param snp_ori Select SNP for the orientation algorithm if orientation=="snp"
#' @param export_order Export the chosen haplotype order after sorting
#' @param import_order Import a given haplotype order
#' @param min_to_plot Only include blocks with at least this many haplotypes in it (default:5)
#' @param intensity Color-intensity for the plotted blocks (default: 0.5 (options: 0-1))
#' @param include Include haplotypes in no block in the orientation position in the graph
#' @param indi number of haplotypes in the dataset
#' @param type length measure (default: "window" , alt: "snp", "bp")
#' @param max_step Maximum distance of blocks to consider in the orientation step (default: 500)
#' @param blocklist block-dataset
#' @param add_sort If FALSE deactivate sorting haplotype for local similarity
#' @param xlim X-axis boundaries of positions to include in the plot
#' @param n_colors Number of different colors used for visualization (without package RColorBrewer limited to 8)
#' @export


plot_block <- function(blocklist, type="snp", orientation="snp", include=TRUE, indi=NULL, min_to_plot = 5,
                      intensity=0.7, add_sort=TRUE, max_step=500,
                      snp_ori=NULL,
                      export_order=FALSE, import_order=FALSE,
                      xlim=NULL, n_colors=20){

  if(length(indi)==0){
    indi <- indi_calc(blocklist)
  }

  orientation_f <- orientation
  se <- blocklist_startend(blocklist, type=type)
  order <- NULL

  if(length(orientation_f)==0){
    order <- 1:indi
  } else if(orientation_f[1]=="front"){
    orientation <- 1:8
  } else if(orientation_f[1]=="mid"){
    orientation <- 1:8 + ceiling(nrow(se)/2-3)
  } else if(orientation_f[1]=="back"){
    orientation <- 1:8 + nrow(se) -8
  } else if(orientation_f[1]=="snp"){
    if(length(snp_ori)==0){
      snp_ori <- mean(c(max(se), min(se)))
    }
    orientation <- which(((se[,1]<= snp_ori)+ (se[,2]>= snp_ori))==2)
  } else{
    order <- orientation_f
  }

  if(length(orientation)>0 && import_order==FALSE){
    sorted <- orientation[1]
    orientation <- orientation[-1]
    while(length(orientation)>0){
      overlap <- numeric(length(orientation))
      actives <- blocklist[[sorted[length(sorted)]]][[6]]
      i1 <- 1
      for(index in orientation){
        nextone <- blocklist[[index]][[6]]
        overlap[i1] <- length(intersect(nextone,actives))/ length(nextone)
        i1 <- i1 + 1
      }
      sorted <- c(sorted, orientation[which.max(overlap)[1]])
      orientation <- orientation[-which.max(overlap)[1]]
    }
    orientation <- sorted

  }
  if(length(order)==0){
    if(include==TRUE){
      orientation <- c(orientation, 0)
    }
    for(index in orientation){
      pl <- length(order)
      if(index==0){
        order <- unique(c(order,1:indi))
        index <- ceiling(median(orientation[-length(orientation)]))
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

      }
    }


  }
  if(length(xlim)>0){
    plot(0,-1000,ylim=c(0,length(order)), ylab="haplotype", xlab="SNP",
         cex.axis=1, cex.lab=1, xlim=xlim)
  } else{
    plot(0,-1000,ylim=c(0,length(order)), xlim=c(1,max(se)), ylab="haplotype", xlab="SNP",
         cex.axis=1, cex.lab=1)
  }


  if(requireNamespace("RColorBrewer")){
    qual_col_pals = RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == 'qual',]
    col_vector = unlist(mapply(RColorBrewer::brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    used_color <- sample(col_vector, n_colors)
  } else{
    used_color <- 1:8
    n_colors <- 8
  }


  activ_colors <- numeric(length(blocklist))
  for(index in 1:length(blocklist)){
    overlap <- duplicated(c(blocklist[[index]][[6]], order))[-(1:blocklist[[index]][[5]])]
    if(sum(overlap) >= min_to_plot){
      taken <- base::intersect(which(se[,1]<=se[index,2]), which(se[,2]>=se[index,1]))
      block_color <- sort(unique(c(0,activ_colors[taken])))
      if(length(block_color)>1 && length(block_color)<=n_colors){
        activ_colors[index] <- sample((1:n_colors)[-block_color],1)
      } else{
        activ_colors[index] <- sample((1:n_colors),1)
      }

      poly <- sort(which(overlap))
      takes <- c(TRUE, diff(poly)>1)
      to_plot <- poly[takes]
      size <- diff(which(takes))
      size <- c(size, length(poly)-sum(size))
      temp1 <- 1

      for(index2 in to_plot){
        polygon(c(se[index,1], se[index,2], se[index,2], se[index,1]), index2-c(1,1,1-size[temp1],1-size[temp1]),
                col=adjustcolor(used_color[activ_colors[index]],alpha.f=intensity), lty=0)
        temp1 <- temp1 +1

      }

    }
  }
  if(export_order==TRUE){
    return(order)
  }
}
