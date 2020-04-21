#' Merge blocks
#'
#' Function to merge/remove/modify blocks
#' @param blocklist block-dataset
#' @param indi number of haplotypes in the dataset
#' @param nwindow number of windows in the dataset
#' @param blockinfo List with all relevant information to each window seperatly
#' @param dhm haploid SNP-dataset
#' @param window_sequence_list sequence of predefined windows (default: NULL ;per row: start$snp, end$snp, length, length - merging_error, start$bp, end$bp)
#' @param dataset dataset which variant nr. for each window
#' @param off_lines minimum number of haplotypes to looose in the creation of a bigger block (default: 5)
#' @param min_similarity minimum rate of the same SNPs to be added to the block (default: 0.99)
#' @param consider_all If TRUE always haplotypes which are not in the node to be in a generated block
#' @param save_allblock If TRUE keep all haplotypes with all windows according to a block (even under min_similarity)
#' @param node_min minimum number of haplotypes per block (default: 5)
#' @param subgroups possible subgroups to consider in the block identification process (default: NULL - list(1:indi))
#' @param min_per_subgroup minimum number of haplotypes per block per subgroup (default: 0)
#' @param subgroup_exception allow for not all subgroups to include a block
#' @param anteil_filter Use only relevant Sequences for filter which haplotypes to check in min_similarity
#' @param intersect_func Used intersect-function (internally relevant for computation time)
#' @param helper Interal information used as a fast abort criterion
#' @param c_dhm Bit-wise coded SNP-dataset
#' @param c_dhm_mode If TRUE use high speed calculation with C in block_merging (default: TRUE)
#' @param run Internal parameter storing the step number
#' @return haplotype library

block_merging <- function(blocklist, blockinfo, dataset, dhm, indi, nwindow, window_sequence_list,  off_lines=5, min_similarity=0.99, consider_all=TRUE,
                          save_allblock=TRUE, node_min=0, subgroups=NULL, min_per_subgroup= 0, anteil_filter=TRUE, helper=NULL,
                          c_dhm=NULL, c_dhm_mode=TRUE, intersect_func=intersect, subgroup_exception=0, run=1){
  if(length(helper)==0){
    helper <- rowSums(blocklist_startend(blocklist, type="snp"))
  } else{
    helper <- rowSums(helper)
  }
  n <- length(helper)

  # Loesche Enthaltende
  be <- blocklist_startend(blocklist, type="snp")
  size <- blocklist_size(blocklist)

  for(index in nrow(be):1){
    possible <- which(((be[index,1]>= be[,1]) * (be[index,2]<= be[,2]) * (1:nrow(be)!=index) * (size>=size[index]))==1)
    delete <- 1
    index2 <- 1
    npos <- length(possible)
    while(index2 < npos && delete==1){
      overlap_blocks <- length(intersect_func(blocklist[[index]][[6]], blocklist[[possible[index2]]][[6]]))
      if( ( overlap_blocks+ off_lines) >=blocklist[[index]][[5]] & overlap_blocks>0){
        delete <- 0
        be[index,1] <- be[index,2] <- 0
        blocklist[[index]] <- "NULL"
      } else{
        index2 <- index2 + 1
      }
    }
  }
  for(index in nrow(be):1){
    if(length(blocklist[[index]])==1){
      blocklist[[index]] <- NULL
    }
  }

  if(run>0){
    if(length(blocklist)>1){
      for(index in length(blocklist):1){
        cluster <- blocklist[[index]][[12]]
        haplotyp <- blocklist[[index]][[4]]

        if(length(blocklist[[index]])>=7 && length(blocklist[[index]][[7]])>0 && length(blocklist[[index]][[7]]$snp)==(blocklist[[index]][[3]]$snp-blocklist[[index]][[2]]$snp+1 )){
          major <- blocklist[[index]][[7]]$snp
        } else {
          major <- numeric(blocklist[[index]][[3]]$snp - blocklist[[index]][[2]]$snp+1)
          for(index2 in blocklist[[index]][[2]]$window:blocklist[[index]][[3]]$window){
            major[window_sequence_list[[cluster]][index2,1]:window_sequence_list[[cluster]][index2,2] - window_sequence_list[[cluster]][blocklist[[index]][[2]]$window,1] +1] <- blockinfo[[cluster]][[index2]][[6]][[haplotyp[index2-blocklist[[index]][[2]]$window+1]]]
          }
          blocklist[[index]][[7]]$snp <- major
        }

      }
    }


    # Kontrolle Enthalten

    for(index in 1:length(blocklist)){
      if(length(blocklist[[index]][[10]])==0 || blocklist[[index]][[10]] != (blocklist[[index]][[3]]$snp - blocklist[[index]][[2]]$snp +1)){

        haplotyp <- blocklist[[index]][[4]]
        major <- blocklist[[index]][[7]]$snp
        anteil <- numeric(indi)
        blockanteil <- numeric(indi)
        cluster <- blocklist[[index]][[12]]

        # Blockanteil-Rechnung
        if(consider_all==TRUE){
          blockanteil <- colSums(dataset[[cluster]][blocklist[[index]][[2]]$window:blocklist[[index]][[3]]$window,, drop=FALSE]==blocklist[[index]][[4]])

        } else{
          blockanteil[blocklist[[index]][[6]]] <- colSums(dataset[[cluster]][blocklist[[index]][[2]]$window:blocklist[[index]][[3]]$window,blocklist[[index]][[6]], drop=FALSE]==blocklist[[index]][[4]])
        }
        # Anteil - Rechnung
        to_consider <- 1:indi
        if(anteil_filter==TRUE){
          to_consider <- to_consider[!(blockanteil==length(haplotyp))]
          if(length(blocklist[[index]])>=9 && length(blocklist[[index]][[9]])>0){
            laenge <- blocklist[[index]][[3]]$snp - blocklist[[index]][[2]]$snp + 1
            max_anteil <- blocklist[[index]][[9]][to_consider,1] + laenge - blocklist[[index]][[9]][to_consider,2]
            to_consider <- to_consider[max_anteil > (laenge * min_similarity)]
          }

        }
        if(!consider_all){
          to_consider <- intersect_func(to_consider, blocklist[[index]][[6]])
        }
        if(length(to_consider)>0){

          if(c_dhm_mode==TRUE){
            anteil[to_consider] <- colSumsEqualSNPs(c_dhm, blocklist[[index]][[2]]$snp, major, to_consider)
          } else{
            anteil[to_consider] <-colSums(dhm[blocklist[[index]][[2]]$snp:blocklist[[index]][[3]]$snp,to_consider, drop=FALSE]==major)

          }
        }

        keep <- which((anteil > ((blocklist[[index]][[3]]$snp - blocklist[[index]][[2]]$snp +1 )*min_similarity)) + save_allblock* (blockanteil == length(haplotyp))>0)


        if(anteil_filter==TRUE){
          if(length(blocklist[[index]])>=9 && length(blocklist[[index]][[9]])>0){
            blocklist[[index]][[9]][to_consider,] <- cbind(anteil[to_consider], blocklist[[index]][[3]]$snp - blocklist[[index]][[2]]$snp + 1)
            blocklist[[index]][[9]][blockanteil==length(haplotyp),] <- blocklist[[index]][[3]]$snp - blocklist[[index]][[2]]$snp + 1
          } else{
            blocklist[[index]][[9]] <- cbind(anteil, blocklist[[index]][[3]]$snp - blocklist[[index]][[2]]$snp +1)
            blocklist[[index]][[9]][blockanteil==length(haplotyp),] <- blocklist[[index]][[3]]$snp - blocklist[[index]][[2]]$snp + 1
          }

        }
        blocklist[[index]][[10]] <- blocklist[[index]][[3]]$snp - blocklist[[index]][[2]]$snp + 1

        previous <- blocklist[[index]][[6]]


        ## Erlaube Statt gleichheit mit Major gleichheit mit einer beliebigen Variante
        #ant <- anteil[blocklist[[index]][[6]]] / (window_size * length(haplotyp))
        #test <- dhm[((blocklist[[index]][[2]]-1)*(window_size)+1):(blocklist[[index]][[3]]*(window_size)),blocklist[[index]][[6]]]==major
        #snp <- dhm[((blocklist[[index]][[2]]-1)*(window_size)+1):(blocklist[[index]][[3]]*(window_size)),blocklist[[index]][[6]]]
        #
        #plot(ksmooth(1:nrow(test), rowSums(test), bandwidth=10))
        ##

        blocklist[[index]][[6]] <- keep

        blocklist[[index]][[5]] <- length(blocklist[[index]][[6]])
      }
    }

  }


  if(min_per_subgroup>0){
    for(index in length(blocklist):1){
      remover <- 0
      exception <- subgroup_exception
      for(group in subgroups){
        if(length(intersect_func(blocklist[[index]][[6]], group))< min_per_subgroup){
          exception <- exception - 1
          if(exception<0){
            remover <- 1
          }
        }
      }
      if(remover==1){
        blocklist[[index]] <- NULL
      }
    }
  }
  blocklist <- blocklist_reorder(blocklist, node_min)



  return(blocklist)
}
