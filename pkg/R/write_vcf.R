#' Function to remove overlapping block segments
#'
#' Function to remove overlapping block segments
#' @param blocklist blocklist
#' @param path path to write the vcf-file to
#' @param chromo Chromosome at hand (default: 1)
#' @param type Position in the vcf file (when available base-bair, alt: "snp", "window")
#' @param hetero Set to TRUE when working on diploid species. Two chromosome sets are assumed to be next to each other in dhm
#' @examples
#' data(ex_maze)
#' blocklist <- block_calculation(ex_maze, bp=1:9999)
#' write_vcf(blocklist, path=tempdir())
#' @export
#' @return VCF-file for the blocklist
#'

write_vcf <- function(blocklist=NULL, path=NULL, chromo = 1, type="bp", hetero = FALSE){

  if(length(path)==0){
    path <- "blocklist.vcf"
  } else{
    path <- paste0(path,".vcf")
  }


  se <- blocklist_startend(blocklist, type=type)
  if(sum(se)==0){
    se <- blocklist_startend(blocklist, type="snp")
  }

  block_dataset <- block_matrix_construction(blocklist)

  if(hetero){
    vcfgeno <- matrix(paste0(block_dataset[,1:(ncol(block_dataset)/2)*1-1], "|", block_dataset[,1:(ncol(block_dataset)/2)*1]), ncol=ncol(block_dataset)/2)
  } else{
    vcfgeno <- matrix(paste0(block_dataset, "|", block_dataset), ncol=ncol(block_dataset))
  }
  ref <- rep("A", nrow(vcfgeno))
  alt <- rep("C", nrow(vcfgeno))

  bp <- se[,1]
  for(index in 2:length(bp)){
    if(bp[index]<= bp[index-1]){
      bp[index] <- bp[index-1]+1
    }
  }

  snpname <- paste0("Block_", 1:nrow(vcfgeno),"(", se[,1],":", se[,2],")")
  options(scipen=999)
  vcfgenofull <- cbind(chromo, as.numeric(bp), snpname, ref, alt, ".", "PASS", ".", "GT", vcfgeno)
  vcfgenofull <- rbind(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", paste0("ID_", 1:ncol(vcfgeno))),vcfgenofull)

  headerfile <- rbind(
    "##fileformat=VCFv4.2",
    gsub("-", "", paste0("##filedate=",  Sys.Date())),
    paste0("##source='HaploBlocker_", utils::packageVersion("HaploBlocker"),"'"),
    "##FORMAT=<ID=GT,Number=1,Type=String,Description='Genotype'>"
  )

  utils::write.table(headerfile, file=path, quote=FALSE, col.names = FALSE, row.names = FALSE)
  utils::write.table(vcfgenofull, file=path, quote=FALSE, col.names = FALSE, row.names = FALSE, append = TRUE, sep="\t")

}
