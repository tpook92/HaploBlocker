#' Calculate Major SNPs
#'
#' Function to calculate the common allele for each SNP of a block
#' @param dhm haploid SNP-dataset
#' @param blocklist blocklist created by blocklist_calculcation algorithm
#' @param recoding If TRUE change allele coding (Major allele "A", Minor allele "C")
#' @return haplotype library with entered allele variants / frequencies

major_snp_calculation <- function(blocklist, dhm, recoding=FALSE){

  if(recoding==TRUE){
    major <- rep("A", nrow(dhm))
    minor <- rep("C", nrow(dhm))
  } else{
    major <- dhm[,1]
    minor <- major
    for(index in 1:nrow(dhm)){
      test <- unique(c(major[index], dhm[index,]))
      if(length(test)>=2){
        minor[index] <- test[2]
      }
    }
  }


  for(index in 1:length(blocklist)){
    snp_seq <- major[blocklist[[index]][[2]]$snp: blocklist[[index]][[3]]$snp]
    freq <- numeric(length(snp_seq))
    count <- rowSums(dhm[blocklist[[index]][[2]]$snp:blocklist[[index]][[3]]$snp, blocklist[[index]][[6]], drop=FALSE]==major[blocklist[[index]][[2]]$snp:blocklist[[index]][[3]]$snp, drop=FALSE])

    freq <- (count/blocklist[[index]][[5]])
    snp_seq[freq<0.5] <- minor[blocklist[[index]][[2]]$snp: blocklist[[index]][[3]]$snp][freq<0.5]

    freq[freq<0.5] <-  1- freq[freq<0.5]
    blocklist[[index]][[7]] <- list()
    blocklist[[index]][[7]]$snp <- snp_seq
    blocklist[[index]][[7]]$freq <- freq
  }



  return(blocklist)
}
