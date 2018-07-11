args <- commandArgs(TRUE)
# args <- c(1,100,"KE", 1)
chromo <- as.numeric(args[1])
chosen_size <- as.numeric(args[2])
line <- args[3]
nsim <- as.numeric(args[4])
target_cov <- as.numeric(args[5])
nr <- as.numeric(args[6])

cov1 <- cov2 <- cov3 <- numeric(nsim)


name <- paste0("RData_MAZE/",line,"_DH_chromo",chromo,".RData")
load(name)


library(HaploBlocker)
library(CHaploBlocker)


for(index in 1:nsim){
  chosen <- sample(1:ncol(data), chosen_size)
  dhm1 <- data[,chosen]
  blocklist <- block_calculation(dhm1, target_coverage = target_cov, max_iteration = 20 , target_stop=0.002)

  cov1[index] <- mean(overlap_test(blocklist, data[,chosen], min_similarity=0.995))
  cov2[index] <- mean(overlap_test(blocklist, data[,-chosen], min_similarity=0.995))
  cov3[index] <- mean(coverage_test(blocklist, chosen_size, type="window"))

  exp_name <- paste0("RData_MAZE/Small_Sample_test/Chr", chromo, "Choosen", chosen_size, "Line", line, "Nsim", nsim, "_tarcov", target_cov,"nr",nr,".RData")
  save(file=exp_name, list=c("cov1", "cov2", "cov3"))
}


