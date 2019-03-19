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
#' @param mindestaehnlichkeit minimum rate of the same SNPs to be added to the block (default: 0.99)
#' @param consider_all If TRUE always haplotypes which are not in the knot to be in a generated block
#' @param save_allblock If TRUE keep all haplotypes with all windows according to a block (even under mindestaehnlichkeit)
#' @param minimum.blocksize minimum.blocksize
#' @param subgroups possible subgroups to consider in the block identification process (default: NULL - list(1:indi))
#' @param min_per_subgroup minimum number of haplotypes per block per subgroup (default: 0)
#' @param anteil_filter Use only relevant Sequences for filter which haplotypes to check in mindestaehnlichkeit
#' @param intersect_func Used intersect-function (internally relevant for computation time)
#' @param helper Interal information used as a fast abort criterion
#' @param c_dhm Bit-wise coded SNP-dataset
#' @param c_dhm_mode If TRUE use high speed calculation with C in block_merging (default: TRUE)
#' @export


block_merging <- function(blocklist, blockinfo, dataset, dhm, indi, nwindow, window_sequence_list,  off_lines=5, mindestaehnlichkeit=0.99, consider_all=TRUE,
                          save_allblock=TRUE, minimum.blocksize=0, subgroups=NULL, min_per_subgroup= 0, anteil_filter=TRUE, helper=NULL,
                          c_dhm=NULL, c_dhm_mode=TRUE, intersect_func=intersect){
  if(length(helper)==0){
    helper <- rowSums(blocklist_startend(blocklist, type="snp"))
  } else{
    helper <- rowSums(helper)
  }
  n <- length(helper)

  if(length(blocklist)>1){
    for(index in length(blocklist):1){
      if(index==1){
        possible_same <- numeric(0)
      } else{
        possible_same <- (1:(index-1))[(helper[1:(index-1)] == helper[index]) * (1:(index-1)) ]
      }

      no.merge <- 1

      cluster <- blocklist[[index]][[12]]
      haplotyp <- blocklist[[index]][[4]]

      if(length(blocklist[[index]])>=7 && length(blocklist[[index]][[7]])>0 && length(blocklist[[index]][[7]]$snp)==(blocklist[[index]][[3]]$snp-blocklist[[index]][[2]]$snp+1 )){
        major <- blocklist[[index]][[7]]$snp
      } else{
        major <- numeric(blocklist[[index]][[3]]$snp - blocklist[[index]][[2]]$snp+1)
        for(index2 in blocklist[[index]][[2]]$window:blocklist[[index]][[3]]$window){
          major[window_sequence_list[[cluster]][index2,1]:window_sequence_list[[cluster]][index2,2] - window_sequence_list[[cluster]][blocklist[[index]][[2]]$window,1] +1] <- blockinfo[[cluster]][[index2]][[6]][[haplotyp[index2-blocklist[[index]][[2]]$window+1]]]
        }
        blocklist[[index]][[7]]$snp <- major
      }

      for(check in possible_same){
        if(no.merge==1 && blocklist[[check]][[2]]$snp == blocklist[[index]][[2]]$snp){

          cluster2 <- blocklist[[check]][[12]]
          haplotyp2 <- blocklist[[check]][[4]]
          if(length(blocklist[[check]])>=7 && length(blocklist[[check]][[7]])>0 && length(blocklist[[check]][[7]]$snp)==(blocklist[[check]][[3]]$snp-blocklist[[check]][[2]]$snp+1 ) ){
            major2 <- blocklist[[check]][[7]]$snp
          } else{
            major2 <- numeric(blocklist[[check]][[3]]$snp - blocklist[[check]][[2]]$snp+1)
            for(index2 in blocklist[[check]][[2]]$window:blocklist[[check]][[3]]$window){
              major2[window_sequence_list[[cluster2]][index2,1]:window_sequence_list[[cluster2]][index2,2] - window_sequence_list[[cluster2]][blocklist[[check]][[2]]$window,1] +1] <- blockinfo[[cluster2]][[index2]][[6]][[haplotyp2[index2-blocklist[[check]][[2]]$window+1]]]
            }
            blocklist[[check]][[7]]$snp <- major2
          }


          merge <- prod(blocklist[[index]][[7]]$snp == blocklist[[check]][[7]]$snp)
          if(blocklist[[check]][[1]][1]=="Off-Window-Variant"){
            merge <- 0
          }
          if(blocklist[[check]][[1]][1]=="Off-Window-Variant" && blocklist[[check]][[1]][1]=="Off-Window-Variant" ){
            if(blocklist[[check]][[2]]$window == blocklist[[index]][[2]]$window && prod( blocklist[[check]][[4]]==blocklist[[index]][[4]])){
              merge <- 1
            }
          }
          if(merge==1 ){
            blocklist[[check]][[6]] <- sort(unique(blocklist[[check]][[6]], blocklist[[index]][[6]]))
            blocklist[[check]][[5]] <- length(blocklist[[check]][[6]])
            blocklist[[index]] <- NULL
            no.merge <-0
          }
        }
      }
    }
  }


  # Loesche Enthaltende
  be <- blocklist_startend(blocklist, type="snp")

  for(index in nrow(be):1){
    possible <- (1:nrow(be))[(be[index,1]>= be[,1]) * (be[index,2]<= be[,2]) *(1:nrow(be)) * (1:nrow(be)!=index)]
    ueberdeckung <- NULL
    for(check in possible){
      control <- prod(blocklist[[check]][[7]]$snp[(be[index,1]-be[check,1]+1):(be[index,2]-be[check,1]+1)] == blocklist[[index]][[7]]$snp)
      if(control==1){
        ueberdeckung <- unique(c(ueberdeckung, blocklist[[check]][[6]]))
      }
    }
    # Sortierung nicht notwendig aber einzige Stelle mit Nicht-Sortierten Intersect-Elementen
    if(length(ueberdeckung)>0){
      ueberdeckung <- sort(ueberdeckung)
    }

    ueberdeckung <- intersect_func(ueberdeckung, blocklist[[index]][[6]])
    if(length(ueberdeckung)>(blocklist[[index]][[5]]-off_lines)){
      blocklist[[index]] <- "NULL"
      be[index,1] <- be[index,2] <- 0
    }
  }
  for(index in nrow(be):1){
    if(length(blocklist[[index]])==1){
      blocklist[[index]] <- NULL
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
          to_consider <- to_consider[max_anteil > (laenge * mindestaehnlichkeit)]
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

      keep <- which((anteil > ((blocklist[[index]][[3]]$snp - blocklist[[index]][[2]]$snp +1 )*mindestaehnlichkeit)) + save_allblock* (blockanteil == length(haplotyp))>0)


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
  if(min_per_subgroup>0){
    for(index in length(blocklist):1){
      remover <- 0
      for(group in subgroups){
        if(length(intersect_func(blocklist[[index]][[6]], group))< min_per_subgroup){
          remover <- 1
        }
      }
      if(remover==1){
        blocklist[[index]] <- NULL
      }
    }
  }
  blocklist <- blocklist_reorder(blocklist, minimum.blocksize)



  return(blocklist)
}
