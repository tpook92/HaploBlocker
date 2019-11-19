#' Optional filtering function
#'
#' Function peform filtering step on the dataset before applying block detection
#' @param dhm haploid SNP-dataset
#' @param maf Minimum minor allel frequency in prefilter (default: 0.05)
#' @param equal_remove If TRUE filter out SNPs in perfect correlation to the next SNP (default: FALSE)
#' @param bp_map X
#' @export

dataset_filter <- function(dhm, maf=0.05, equal_remove=FALSE, bp_map=NULL){
  maf1 <- numeric(nrow(dhm))
  for(index in 1:length(maf1)){
    maf1[index] <- sum(dhm[index,]==dhm[index,1])
  }
  maf1 <- maf1 / ncol(dhm)
  maf1[maf1>0.5] <- 1-maf1[maf1>0.5]
  removes <- which(maf1<maf)

  if(length(removes)>0){
    dhm <- dhm[-removes,]
    maf1 <- maf1[-removes]
    bp_map <- bp_map[-removes]
  }
  removes2 <- NULL
  if(equal_remove==TRUE){
    for(index in (length(maf1)):2){
      if(maf1[index]==maf1[index-1]){
        first <- which(dhm[index,]==dhm[index,1])
        second <- which(dhm[index-1,]==dhm[index-1,1])
        if(length(second)!=length(first)){
          second <- which(dhm[index-1,]!=dhm[index-1,])
        }
        if(sum(first==second)==length(first)){
          removes2 <- c(removes2, index)
        }

      }
    }
    dhm <- dhm[-removes2,]
    maf1 <- maf1[-removes2]
    bp_map <- bp_map[-removes2]
  }
  return(list(dhm, bp_map))
}
