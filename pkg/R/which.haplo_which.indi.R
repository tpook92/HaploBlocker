#' Which block
#'
#' Function to derivde in which block a giving haplotype is
#' @param blocklist block-dataset
#' @param nr number of haplotype to analyse
#' @examples
#' data(blocklist_ex_maze)
#' which.block(blocklist_ex_maze, 1)
#' @export
#' @return Binary vector coding present/absense of a haplotype in each block

which.block <- function(blocklist, nr){
  include <- numeric(length(blocklist))
  for(index in 1:length(blocklist)){
    if(sum(blocklist[[index]][[6]]==nr)>0){
      include[index] <-1
    }

  }
  return(include)
}

#' Which indi
#'
#' Function to derive in which haplotypes are locally in the same block as a giving haplotype
#' @param blocklist block-dataset
#' @param nr number of haplotype to analyse
#' @param start Lower boundary of the considered intervall
#' @param end Upper boundary of the considered intervall
#' @param se start/end points of blocks (blocklist_start_end() output)
#' @param type length measure (default: "bp" , alt: "snp", "window")
#' @examples
#' data(blocklist_ex_maze)
#' which.indi(blocklist_ex_maze, 1, 1,1, type="window")
#' @export
#' @return List of haplotype that share same local haplotype blocks

which.indi <- function(blocklist, nr, start, end, se=NULL, type="bp"){
  indis <- NULL
  if(length(se)==0){
    se <- blocklist_startend(blocklist, type=type)
  }

  for(index in (1:length(blocklist))[se[,1]<= start & se[,2]>= end]){
    if(sum(blocklist[[index]][[6]]==nr)>0){
      indis <- sort(unique(c(indis, blocklist[[index]][[6]])))
    }
  }
  return(indis)
}
