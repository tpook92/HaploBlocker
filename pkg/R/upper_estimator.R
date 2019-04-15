#' Library presents diploid
#'
#' Calculate upper limited of present in a blocklist (ignore possible phasing errors)
#' @param blocklist block-dataset
#' @param s0m dataset of diploid individuals to check
#' @param min_similarity XXX
#' @export

diploid_upper <- function(blocklist, s0m, min_similarity=0.99){
  n_block <- length(blocklist)
  max_size <- numeric(n_block)

  for(index in 1:n_block){
    start <- blocklist[[index]][[2]]$snp
    end <- blocklist[[index]][[3]]$snp
    seq <- blocklist[[index]][[7]]$snp
    same <- s0m[start:end,]==seq
    for(index2 in 1:(ncol(s0m)/2)){
      same[,index2] <- (same[,index2*2] + same[,index2*2-1])>0
    }
    max_sim <- colSums(same[,1:index2]) / nrow(same)
    max_size[index] <- sum(max_sim>=min_similarity)

  }
  return(max_size)
}

#' Library presents diploid
#'
#' Calculate upper limited of present in a blocklist (ignore possible phasing errors)
#' @param blocklist block-dataset
#' @param s0m dataset of diploid individuals to check
#' @param min_similarity XXX
#' @export

diploid_upper2 <- function(blocklist, s0m, min_similarity=0.99){
  n_block <- length(blocklist)
  max_size <- numeric(n_block)
  max_double <- numeric(n_block)
  for(index in 1:n_block){
    start <- blocklist[[index]][[2]]$snp
    end <- blocklist[[index]][[3]]$snp
    seq <- blocklist[[index]][[7]]$snp
    same <- s0m[start:end,]==seq
    double_check <- colSums(same) / nrow(same)
    max_double <- numeric((ncol(s0m)/2))
    for(index2 in 1:(ncol(s0m)/2)){
      same[,index2] <- (same[,index2*2] + same[,index2*2-1])>0
      max_double[index2] <- sum(double_check[c(index2*2, index2*2-1)]>=min_similarity)
    }
    max_sim <- colSums(same[,1:(ncol(s0m)/2)]) / nrow(same)
    maxt <- max_single <- (max_sim>min_similarity)
    maxt[max_double>max_single] <- max_double[max_double>max_single]
    max_size[index] <- sum(maxt)

  }
  return(max_size)
}

#' Library presents haploid
#'
#' Calculate upper limited of present in a blocklist (not-ignore possible phasing errors)
#' @param blocklist block-dataset
#' @param s0m dataset of diploid individuals to check
#' @param min_similarity XXX
#' @export

haploid_upper <- function(blocklist, s0m, min_similarity=0.99){
  n_block <- length(blocklist)
  max_size <- numeric(n_block)

  for(index in 1:n_block){
    start <- blocklist[[index]][[2]]$snp
    end <- blocklist[[index]][[3]]$snp
    seq <- blocklist[[index]][[7]]$snp
    same <- s0m[start:end,]==seq
    max_sim <- colSums(same) / nrow(same)
    max_size[index] <- sum(max_sim>=min_similarity)

  }
  return(max_size)

}
