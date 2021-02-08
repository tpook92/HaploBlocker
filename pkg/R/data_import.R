#' Function to import genomic data from a vcf/ped file
#'
#' Function to import genomic data from a vcf/ped file
#' @param dhm Path of the file to import
#' @param inbred Set to TRUE when working with inbred material and importing genomic data from file
#' @return genomic dataset and basepair positions from the vcf/ped-file
#' @export
#'

data_import <- function(dhm, inbred=FALSE, verbose=TRUE){

  data_type <- substr(dhm, start= nchar(dhm)-3, stop= nchar(dhm))

  if(data_type==".vcf" || "f.gz"){
    if(verbose) cat("Data input identified as vcf-file - extract genomic information. \n")
    if(requireNamespace("vcfR", quietly = TRUE)){
      vcf_file <- vcfR::read.vcfR(dhm)
      haplo1 <- substr(vcf_file@gt[,-1], start=1, stop=1)
      haplo2 <- substr(vcf_file@gt[,-1], start=3, stop=3)
      haplo <- cbind(haplo1, haplo2)
      bp <- as.numeric(vcf_file@fix[,2])
    } else{
      vcf_file <- utils::read.table(dhm)
      vcf_data <- as.matrix(vcf_file[,-(1:9)])
      haplo1 <- substr(vcf_data, start=1, stop=1)
      haplo2 <- substr(vcf_data, start=3, stop=3)
      bp <- as.numeric(vcf_file[,2])
    }
    if(mean(haplo1==haplo2)>0.95){
      cat("Your material seems to be highly homozygous. Consider setting the parameter inbred to TRUE (only consider 1 haplotype per line)!\n")
      warnings("Your material seems to be highly homozygous. Consider setting the parameter inbred to TRUE (only consider 1 haplotype per line)!")
    }

    if(inbred){
      haplo <- haplo1
    } else{
      haplo <- haplo[,c(0,ncol(haplo1)) + rep(1:ncol(haplo1), each=2)]
    }



  } else if(data_type==".ped" || data_type=="d.gz"){
    bp <- NULL
    if(verbose) cat("Data input identified as Ped-map-file - extract genomic information. \n")
    if(verbose) cat("Haplotype phase is assumed to be by colum - No internal phasing performed! \n")
    ped_file <- utils::read.table(dhm)
    haplo12 <- t(ped_file[,-(1:6)])
    haplo <- matrix(0, ncol = ncol(haplo12)*2, nrow=nrow(haplo12)/2)
    for(index1 in 1:ncol(haplo12)){
      haplo[,index1*2+c(-1,0)] <- matrix(haplo12[,index1], ncol=2, byrow=TRUE)
    }

    if(inbred){
      haplo <- haplo[,(1:(ncol(haplo)/2))*2]
    }
    if( mean(haplo[,(1:(ncol(haplo)/2))*2] == haplo[,(1:(ncol(haplo)/2))*2-1])>0.95){
      cat("Your material seems to be highly homozygous. Consider setting the parameter inbred to TRUE (only consider 1 haplotype per line)!\n")
      warnings("Your material seems to be highly homozygous. Consider setting the parameter inbred to TRUE (only consider 1 haplotype per line)!")

    }
  } else{
    stop("Data type could not be identified. Please manually import the dataset. \n")
  }
  #storage.mode(haplo) <- "integer"

  return(list(haplo, bp))

}
