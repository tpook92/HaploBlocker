#' Construct blockinfo dataset
#'
#' Internal Function to generate a dataset containing which variant in present in with window
#' @param blockinfo List with all relevant information to each window seperatly
#' @param blocklist blocklist created by blocklist_calculcation algorithm
#' @param indi number of haplotypes in the dataset
#' @param nwindow number of windows in the dataset
#' @export

block_dataset_construction <- function(blockinfo, blocklist=NULL, indi=NULL, nwindow=NULL){
  if(length(indi)==0){
    indi <- indi_calc(blocklist)
  }
  if(length(nwindow)==0){
    nwindow <- length(blockinfo)
  }
  dataset <- matrix(0, nrow=nwindow, ncol=indi)
  for(index in 1:length(blockinfo)){
    for(index2 in 1:length(blockinfo[[index]][[5]])){
      dataset[index, blockinfo[[index]][[5]][[index2]]] <- index2
    }
  }
  return(dataset)
}

#' Construct blocklist dataset
#'
#' Internal Function to generate a dataset containing which variant in present in with window
#' @param blocklist blocklist created by blocklist_calculcation algorithm
#' @param indi number of haplotypes in the dataset
#' @export

block_matrix_construction <- function(blocklist, indi=NULL){
  if(length(indi)==0){
    indi <- indi_calc(blocklist)
  }
  dataset <- matrix(0, nrow=length(blocklist), ncol=indi)
  for(index in 1:length(blocklist)){
    dataset[index,blocklist[[index]][[6]]] <- 1
    }

  return(dataset)
}
