#' Function to derive bEHH
#'
#' Function to generate haplotype blocks from haploid data
#' @param blocklist blocklist
#' @param data haploid SNP-dataset
#' @param marker marker for which to calculate eHH/iHH
#' @param plot generate a bEHH-Plot (default: FALSE)
#' @param position1 Position in bp (default: 1,2,3,...)
#' @param standardization Standardization by allele freq (0), allele freq of individuals in blocks (1), block (2), non (3)
#' @param group List containing all individuals assign in each group to derive different scores for subgroups
#' @param return_ehh Return single ehh values
#' @param verbose Set to FALSE to not display any prints
#' @examples
#' data(blocklist_ex_maze)
#' ehh_scores <- block_ehh(blocklist_ex_maze, marker = 5000, plot=TRUE)
#' @export
#' @return Block-based EHH scores


block_ehh <- function(blocklist=NULL, data=NULL, marker, plot=FALSE, position1=NULL, standardization=3, group=NULL,
                      return_ehh=TRUE, verbose = TRUE){


  if(length(blocklist)==0){
    blocklist <- block_calculation(data, verbose = verbose)
  }
  if(length(position1)==0){
    position1 <- 1:max(blocklist_startend(blocklist))
  }
  if(length(group)==0){
    group <- list()
    group[[1]] <- 1:indi_calc(blocklist)
  }

  iHH_list <- NULL
  ehh_list <- NULL

  for(grp in 1:length(group)){
    blockl <- blocklist
    se <- blocklist_startend(blockl)
    active_blocks <- which(((se[,1]<=marker) * (se[,2]>=marker))==1)
    if(length(active_blocks)==0){
      return(0)
    }
    marker_low <- max(se[active_blocks,1])
    marker_high <- min(se[active_blocks,2])
    window <- c(min(se[active_blocks,1]), max(se[active_blocks,2]))
    variants <- matrix(NA, nrow=length(active_blocks), ncol=diff(window)+1)

    sub <- window[1] -1
    for(index in 1:length(active_blocks)){
      start <- blockl[[active_blocks[[index]]]][[2]]$snp
      end <- blockl[[active_blocks[[index]]]][[3]]$snp
      variants[index, (start:end)-sub] <- blockl[[active_blocks[[index]]]][[7]]$snp
    }
    variants[is.na(variants)] <- 256 # max 255 variants

    boundaries <- sort(unique(c(marker_low, (se[active_blocks,2]))))[-1]
    variant <- list()
    class <- list()
    class_count <- list()
    for(z in 1:length(boundaries)){
      variant[[z]] <- unique(variants[,marker_low:boundaries[z]-sub, drop=FALSE])
      variant[[z]] <- variant[[z]][variant[[z]][,boundaries[z]-marker_low+1]!=256,,drop=FALSE]
      class[[z]] <- list()
      class_count[[z]] <- rep(0, nrow(variant[[z]]))
      for(index in 1:nrow(variants)){
        assign_g <- which.max(colSums(t(variant[[z]])==variants[index,marker_low:boundaries[z]-sub]))
        if(variants[index,boundaries[z]-sub]!=256){
          if(length(class[[z]])<assign_g){
            class[[z]][[assign_g]] <- blockl[[active_blocks[index]]][[6]]
          } else{
            class[[z]][[assign_g]] <- unique(c(class[[z]][[assign_g]], blockl[[active_blocks[index]]][[6]]))
          }
        }


      }

      for(index in 1:nrow(variant[[z]])){
        if(length(class[[z]])>=index){
          class_count[[z]][index] <- length( base::intersect(class[[z]][[index]], group[[grp]]))
        }

      }

      while(sum(duplicated(unlist(class[[z]]))>0)){

        full <- unlist(class[[z]])
        full <- full[which(duplicated(full))]
        candidates_remove <- NULL
        for(index in 1:length(class[[z]])){
          if(sum(class[[z]][[index]]==full[1])){
            candidates_remove <- c(candidates_remove, index)
          }
        }
        candidates_remove <- candidates_remove[-which.max(class_count[[z]][candidates_remove])]
        for(index in candidates_remove){
          class[[z]][[index]] <- setdiff(class[[z]][[index]], full)
        }
        for(index in 1:nrow(variant[[z]])){
          class_count[[z]][index] <- length( base::intersect(class[[z]][[index]], group[[grp]]))
        }

      }



    }

    if(standardization==0){
      allel <- c(sum(data[marker,group[[grp]]]==0), sum(data[marker,group[[grp]]]==2))
      ehh_divisor <- sum(choose(allel,2))
    } else if(standardization==1){
      allel <- class_count[[1]]
      store <- allel
      ehh_divisor <- sum(choose(allel,2))
    } else if(standardization==2){
      allel <- class_count[[2]]
      ehh_divisor <- sum(choose(allel,2))
    } else{
      ehh_divisor <- choose(length(group[[grp]]),2)
    }

    ehh_dividend <- numeric(length(class_count))
    for(index in 1:length(ehh_dividend)){
      ehh_dividend[index] <- sum(choose(class_count[[index]],2))
    }
    value <- ehh_dividend/ehh_divisor



    position <- position1[boundaries]


    boundaries2 <- sort(unique(c(marker_high, (se[active_blocks,1]))), decreasing = TRUE)[-1]
    variant <- list()
    class <- list()
    class_count <- list()
    for(z in 1:length(boundaries2)){
      variant[[z]] <- unique(variants[,boundaries2[z]:marker_high-sub, drop=FALSE])
      variant[[z]] <- variant[[z]][variant[[z]][,1]!=256,,drop=FALSE]
      class[[z]] <- list()
      class_count[[z]] <- rep(0, nrow(variant[[z]]))
      for(index in 1:nrow(variants)){
        assign_g <- which.max(colSums(t(variant[[z]])==variants[index,boundaries2[z]:marker_high-sub]))
        if(variants[index,boundaries2[z]-sub]!=256){
          if(length(class[[z]])<assign_g){
            class[[z]][[assign_g]] <- blockl[[active_blocks[index]]][[6]]
          } else{
            class[[z]][[assign_g]] <- unique(c(class[[z]][[assign_g]], blockl[[active_blocks[index]]][[6]]))
          }
        }


      }
      for(index in 1:nrow(variant[[z]])){
        if(length(class[[z]])>=index){
          class_count[[z]][index] <- length( base::intersect(class[[z]][[index]], group[[grp]]))
        }

      }
      while(sum(duplicated(unlist(class[[z]]))>0)){

        full <- unlist(class[[z]])
        full <- full[which(duplicated(full))]
        candidates_remove <- NULL
        for(index in 1:length(class[[z]])){
          if(sum(class[[z]][[index]]==full[1])){
            candidates_remove <- c(candidates_remove, index)
          }
        }
        candidates_remove <- candidates_remove[-which.max(class_count[[z]][candidates_remove])]
        for(index in candidates_remove){
          class[[z]][[index]] <- setdiff(class[[z]][[index]], full)
        }
        for(index in 1:nrow(variant[[z]])){
          class_count[[z]][index] <- length( base::intersect(class[[z]][[index]], group[[grp]]))
        }

      }


    }

    if(standardization==0){
      allel <- c(sum(data[marker,group[[grp]]]==0), sum(data[marker,group[[grp]]]==2))
      ehh_divisor <- sum(choose(allel,2))
    } else if(standardization==1){
      allel <- store
      ehh_divisor <- sum(choose(allel,2))
    } else if(standardization==2){
      allel <- class_count[[1]]
      ehh_divisor <- sum(choose(allel,2))
    } else{
      ehh_divisor <- choose(length(group[[grp]]),2)
    }
    ehh_dividend <- numeric(length(class_count))
    for(index in 1:length(ehh_dividend)){
      ehh_dividend[index] <- sum(choose(class_count[[index]],2))
    }

    if(length(ehh_dividend)>1){
      value2 <- c((ehh_dividend/ehh_divisor)[length(ehh_dividend):2], value)
    } else{
      value2 <- value
    }
    value <- c((ehh_dividend/ehh_divisor)[length(ehh_dividend):1], value)


    position <-c(position1[boundaries2][length(boundaries2):1], position)

    if(standardization!=0 & standardization!=1){
      value2[which(diff(position)==0)] <- value2[which(diff(position)==0)+1]
    }

    if(plot){
      plot(0,-10, ylim=c(0,max(1, value2)), xlim=c(min(position), max(position)), xlab="Position (Mb)", ylab="Block EHH")

      for(index in 1:length(position)){
        lines(position[index:(index+1)], rep(value2[index],2))
      }
      for(index in 1:length(position)){
        lines(rep(position[(index+1)],2), value2[index:(index+1)])
      }
      abline(v=marker, col="red", lwd=2)
    }
    iHS <- sum(value2 * diff(position))

    iHH_list <- c(iHH_list, iHS)
    if(return_ehh){
      ehh <- numeric(length(position1))
      for(index in 1:(length(position)-1)){
        ehh[position[index]:position[index+1]] <- value2[index]
      }
      ehh_list <- rbind(ehh_list, ehh)
    }
  }

  if(return_ehh){
    return(ehh_list)
  } else{
    return(iHH_list)
  }

}
