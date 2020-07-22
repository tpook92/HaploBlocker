#' Mainfunction to calculate haplotype blocks
#'
#' Function to generate haplotype blocks from haploid data
#' @param blocklist blocklist
#' @param data haploid SNP-dataset
#' @param plot generate a IHH-Plot (default: FALSE)
#' @param position1 Position in bp (default: 1,2,3,...)
#' @param standardization Standardization by allele freq (0), allele freq of individuals in blocks (1), block (2), non (3)
#' @param group List containing all individuals assign in each group to derive different scores for subgroups
#' @param verbose Set to FALSE to not display any prints
#' @examples
#' data(blocklist_ex_maze)
#' ihh_scores <- block_ihh(blocklist_ex_maze, plot=TRUE)
#' @export
#' @return Block-based IHH scores

block_ihh <- function(blocklist=NULL, data=NULL, plot=FALSE, position1=NULL, standardization=3,
                      group=NULL, verbose = TRUE){

  if(length(group)==0){
    group <- list()
    group[[1]] <- 1:indi_calc(blocklist)
  }
  if(length(blocklist)==0){
    library(HaploBlocker)
    blocklist <- block_calculation(data, verbose = verbose)
  }
  if(length(position1)==0){
    position1 <- 1:max(blocklist_startend(blocklist))
  }

  se <- blocklist_startend(blocklist, type="snp")

  end_block  <- sort(unique(c(0,se[,1]-1, se[,2])))[-1]
  start_block <- c(min(se), end_block[1:(length(end_block)-1)]+1)

  ihh_scores <- matrix(NA, nrow=length(group), ncol=length(position1))
  for(index in 1:length(start_block)){
    ihh_scores[,start_block[index]:end_block[index]] <-  block_ehh(blocklist, data=data, marker = start_block[index], position1 = position1,
                                                                   standardization = standardization, return_ehh=FALSE, group=group)
  }

  ihh_scores[is.na(ihh_scores)] <- 0
  if(plot){

    plot(position1,ihh_scores[1,], type="l", ylab="IHH score", xlab="position", ylim=c(0, max(ihh_scores)))
    if(nrow(ihh_scores)>1){
      for(index in 2:nrow(ihh_scores)){
        lines(position1,ihh_scores[index,], col=index)
      }


    }
  }

  return(ihh_scores)

}
