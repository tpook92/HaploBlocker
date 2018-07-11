
name <- paste0("C:/Users/tpook/Desktop/Genetic_Datasets/Batch3_KEPE/","KE","_DH_chromo",1,".RData")


load(name)
library(HaploBlocker)
data1 <- data[1:20000,]
blocklist <- block_calculation(data1)


png(file="C:/Users/tpook/Desktop/Fig4.png", width=2250, height=1650, res=300)
par(ps = 12, cex = 1, cex.main = 1, cex.lab=1, cex.axis=1)
plot_block(blocklist, indi=501)
dev.off()

{
  png(file="C:/Users/tpook/Desktop/Fig5.png", width=2250, height= 1360, res=300)
  race <- "KE"
  diff <- 1
  closer <- -0.5
  par(ps = 12, cex = 1, cex.main = 1, cex.lab=1, cex.axis=1)
  par(mfrow=c(1,2))
  par(mar=c(4.6,3.1,2.6,0.1))
  zwei <- c(4.6,2.6,2.6,0.6)
#  par(mar=c(4.6,4.1,2.6,0.1))
#  zwei <- c(4.6,2.6,2.6,1.6)
  scov1m <- scov2m <- scov3m <- numeric(8)
  up1 <- up2 <- up3 <- down1 <- down2 <- down3 <- numeric(8)
  for(chromo in 1:1){
    sample_range <- c(30,50,75,100,150,200,300,400)
    scov1 <- NULL
    scov2 <- NULL
    scov3 <- NULL
    sd1 <- NULL
    sd2 <- NULL
    sd3 <- NULL

    for(sample in sample_range){
      name <- paste0("C:/Users/tpook/Desktop/Small_Sample_test/Chr", chromo, "Choosen",sample, "Line",race,"Nsim10.RData")
      load(name)
      scov1 <- cbind(scov1, cov1t)
      scov2 <- cbind(scov2, cov2t)
      scov3 <- cbind(scov3, cov3t)
    }
  }
  scov1[scov1==0] <- NA
  scov2[scov2==0] <- NA
  scov3[scov3==0] <- NA

  scov4 <- cbind(scov2[,1], scov3[,1], scov2[,2],scov3[,2],
                 scov2[,3], scov3[,3], scov2[,4],scov3[,4],
                 scov2[,5], scov3[,5], scov2[,6],scov3[,6],
                 scov2[,7], scov3[,7], scov2[,8],scov3[,8])


  ats <- sort(rep(1:8,2))
  boxplot(scov4, at=ats,
          col=rep(c(2,8),8), outcol=rep(c(2,1),8), xaxt="n", ylab="",xlab="training set (of 501)",
          ylim=c(0.43,0.94), main="Default settings:")
  title(ylab="coverage", line=2.2)
  axis(side=1, at=c(1+seq(0,by=diff,length.out = 8))[1:4*2-1], labels=sample_range[1:4*2-1],
       xlab="training set (of 501)")
  axis(side=1, at=c(1+seq(0,by=diff,length.out = 8))[1:4*2], labels=sample_range[1:4*2],
       xlab="training set (of 501)")
  par(mar=zwei)
  scov1m <- scov2m <- scov3m <- numeric(8)
  up1 <- up2 <- up3 <- down1 <- down2 <- down3 <- numeric(8)
  for(chromo in 1:1){
    sample_range <- c(30,50,75,100,150,200,300,400)
    scov1 <- NULL
    scov2 <- NULL
    scov3 <- NULL
    sd1 <- NULL
    sd2 <- NULL
    sd3 <- NULL

    for(sample in sample_range){
      name <- paste0("C:/Users/tpook/Desktop/Small_Sample_test/Chr", chromo, "Choosen",sample, "Line",race,"Nsim10_tarcov0.9.RData")
      load(name)
      scov1 <- cbind(scov1, cov1t)
      scov2 <- cbind(scov2, cov2t)
      scov3 <- cbind(scov3, cov3t)
    }
  }
  scov1[scov1==0] <- NA
  scov2[scov2==0] <- NA
  scov3[scov3==0] <- NA

  scov4 <- cbind(scov2[,1], scov3[,1], scov2[,2],scov3[,2],
                 scov2[,3], scov3[,3], scov2[,4],scov3[,4],
                 scov2[,5], scov3[,5], scov2[,6],scov3[,6],
                 scov2[,7], scov3[,7], scov2[,8],scov3[,8])


  ats <- sort(rep(1:8,2))
  boxplot(scov4, at=ats,
          col=rep(c(2,8),8), outcol=rep(c(2,1),8), xaxt="n", ylab="",xlab="training set (of 501)",
          ylim=c(0.43,0.94), main="Fixed coverage in training set:")
  axis(side=1, at=c(1+seq(0,by=diff,length.out = 8))[1:4*2-1], labels=sample_range[1:4*2-1],
       xlab="training set (of 501)")
  axis(side=1, at=c(1+seq(0,by=diff,length.out = 8))[1:4*2], labels=sample_range[1:4*2])
  legend("bottomright", c("training set", "test set"), fill=c(8,2))

  dev.off()
}

load("C:/Users/tpook/Desktop/R-Stuff/b123.RData")

png(file="C:/Users/tpook/Desktop/Fig6.png", width=1559, height=1245, res=300)

par(ps = 12, cex = 1, cex.main = 1, cex.lab=1, cex.axis=1)
par(mar=c(5.1,4.1,2.1,2.1))
ct <- c(0.51,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1)
type1 <- type2 <- type1b <- type2b <-  ci1 <- ci2 <- numeric(length(ct))

type1[1] <- 1 - (sum(b051_total[,3]) / sum(b051_total[,1]))
type2[1] <- (sum(b051_total[,2])-sum(b051_total[,3])) / (50*160*501 - sum(b051_total[,4]))
ci1[1] <- sd(b051_total[,3]/b051_total[,1])/sqrt(nrow(b051_total)) * 1.96
ci2[1] <- sd((b051_total[,2]-b051_total[,3])/(501-b051_total[,4]))/sqrt(nrow(b051_total)) * 1.96
type1b[1] <- 1 - (sum(b051_total[,6]) / sum(b051_total[,4]))
type2b[1] <- (sum(b051_total[,5])-sum(b051_total[,6])) / (50*160*501 - sum(b051_total[,4]))

type1[2] <- 1 - (sum(b055_total[,3]) / sum(b055_total[,1]))
type2[2] <- (sum(b055_total[,2])-sum(b055_total[,3])) / (50*160*501 - sum(b055_total[,4]))
ci1[2] <- sd(b055_total[,3]/b055_total[,1])/sqrt(nrow(b055_total)) * 1.96
ci2[2] <- sd((b055_total[,2]-b055_total[,3])/(501-b055_total[,4]))/sqrt(nrow(b055_total)) * 1.96
type1b[2] <- 1 - (sum(b055_total[,6]) / sum(b055_total[,4]))
type2b[2] <- (sum(b055_total[,5])-sum(b055_total[,6])) / (50*160*501 - sum(b055_total[,4]))

type1[3] <- 1 - (sum(b060_total[,3]) / sum(b060_total[,1]))
type2[3] <- (sum(b060_total[,2])-sum(b060_total[,3])) / (50*160*501 - sum(b060_total[,4]))
ci1[3] <- sd(b060_total[,3]/b060_total[,1])/sqrt(nrow(b060_total)) * 1.96
ci2[3] <- sd((b060_total[,2]-b060_total[,3])/(501-b060_total[,4]))/sqrt(nrow(b060_total)) * 1.96
type1b[3] <- 1 - (sum(b060_total[,6]) / sum(b060_total[,4]))
type2b[3] <- (sum(b060_total[,5])-sum(b060_total[,6])) / (50*160*501 - sum(b060_total[,4]))

type1[4] <- 1 - (sum(b065_total[,3]) / sum(b065_total[,1]))
type2[4] <- (sum(b065_total[,2])-sum(b065_total[,3])) / (50*160*501 - sum(b065_total[,4]))
ci1[4] <- sd(b065_total[,3]/b065_total[,1])/sqrt(nrow(b065_total)) * 1.96
ci2[4] <- sd((b065_total[,2]-b065_total[,3])/(501-b065_total[,4]))/sqrt(nrow(b065_total)) * 1.96
type1b[4] <- 1 - (sum(b065_total[,6]) / sum(b065_total[,4]))
type2b[4] <- (sum(b065_total[,5])-sum(b065_total[,6])) / (50*160*501 - sum(b065_total[,4]))

type1[5] <- 1 - (sum(b070_total[,3]) / sum(b070_total[,1]))
type2[5] <- (sum(b070_total[,2])-sum(b070_total[,3])) / (50*160*501 - sum(b070_total[,4]))
ci1[5] <- sd(b070_total[,3]/b070_total[,1])/sqrt(nrow(b070_total)) * 1.96
ci2[5] <- sd((b070_total[,2]-b070_total[,3])/(501-b070_total[,4]))/sqrt(nrow(b070_total)) * 1.96
type1b[5] <- 1 - (sum(b070_total[,6]) / sum(b070_total[,4]))
type2b[5] <- (sum(b070_total[,5])-sum(b070_total[,6])) / (50*160*501 - sum(b070_total[,4]))

type1[6] <- 1 - (sum(b075_total[,3]) / sum(b075_total[,1]))
type2[6] <- (sum(b075_total[,2])-sum(b075_total[,3])) / (50*160*501 - sum(b075_total[,4]))
ci1[6] <- sd(b075_total[,3]/b075_total[,1])/sqrt(nrow(b075_total)) * 1.96
ci2[6] <- sd((b075_total[,2]-b075_total[,3])/(501-b075_total[,4]))/sqrt(nrow(b075_total)) * 1.96
type1b[6] <- 1 - (sum(b075_total[,6]) / sum(b075_total[,4]))
type2b[6] <- (sum(b075_total[,5])-sum(b075_total[,6])) / (50*160*501 - sum(b075_total[,4]))

type1[7] <- 1 - (sum(b080_total[,3]) / sum(b080_total[,1]))
type2[7] <- (sum(b080_total[,2])-sum(b080_total[,3])) / (50*160*501 - sum(b080_total[,4]))
ci1[7] <- sd(b080_total[,3]/b080_total[,1])/sqrt(nrow(b080_total)) * 1.96
ci2[7] <- sd((b080_total[,2]-b080_total[,3])/(501-b080_total[,4]))/sqrt(nrow(b080_total)) * 1.96
type1b[7] <- 1 - (sum(b080_total[,6]) / sum(b080_total[,4]))
type2b[7] <- (sum(b080_total[,5])-sum(b080_total[,6])) / (50*160*501 - sum(b080_total[,4]))

type1[8] <- 1 - (sum(b085_total[,3]) / sum(b085_total[,1]))
type2[8] <- (sum(b085_total[,2])-sum(b085_total[,3])) / (50*160*501 - sum(b085_total[,4]))
ci1[8] <- sd(b085_total[,3]/b085_total[,1])/sqrt(nrow(b085_total)) * 1.96
ci2[8] <- sd((b085_total[,2]-b085_total[,3])/(501-b085_total[,4]))/sqrt(nrow(b085_total)) * 1.96
type1b[8] <- 1 - (sum(b085_total[,6]) / sum(b085_total[,4]))
type2b[8] <- (sum(b085_total[,5])-sum(b085_total[,6])) / (50*160*501 - sum(b085_total[,4]))

type1[9] <- 1 - (sum(b090_total[,3]) / sum(b090_total[,1]))
type2[9] <- (sum(b090_total[,2])-sum(b090_total[,3])) / (50*160*501 - sum(b090_total[,4]))
ci1[9] <- sd(b090_total[,3]/b090_total[,1])/sqrt(nrow(b090_total)) * 1.96
ci2[9] <- sd((b090_total[,2]-b090_total[,3])/(501-b090_total[,4]))/sqrt(nrow(b090_total)) * 1.96
type1b[9] <- 1 - (sum(b090_total[,6]) / sum(b090_total[,4]))
type2b[9] <- (sum(b090_total[,5])-sum(b090_total[,6])) / (50*160*501 - sum(b090_total[,4]))

type1[10] <- 1 - (sum(b095_total[,3]) / sum(b095_total[,1]))
type2[10] <- (sum(b095_total[,2])-sum(b095_total[,3])) / (50*160*501 - sum(b095_total[,4]))
ci1[10] <- sd(b095_total[,3]/b095_total[,1])/sqrt(nrow(b095_total)) * 1.96
ci2[10] <- sd((b095_total[,2]-b095_total[,3])/(501-b095_total[,4]))/sqrt(nrow(b095_total)) * 1.96
type1b[10] <- 1 - (sum(b095_total[,6]) / sum(b095_total[,4]))
type2b[10] <- (sum(b095_total[,5])-sum(b095_total[,6])) / (50*160*501 - sum(b095_total[,4]))

type1[11] <- 1 - (sum(b100_total[,3]) / sum(b100_total[,1]))
type2[11] <- (sum(b100_total[,2])-sum(b100_total[,3])) / (50*160*501 - sum(b100_total[,4]))
ci1[11] <- sd(b100_total[,3]/b100_total[,1])/sqrt(nrow(b100_total)) * 1.96
ci2[11] <- sd((b100_total[,2]-b100_total[,3])/(501-b100_total[,4]))/sqrt(nrow(b100_total)) * 1.96
type1b[11] <- 1 - (sum(b100_total[,6]) / sum(b100_total[,4]))
type2b[11] <- (sum(b100_total[,5])-sum(b100_total[,6])) / (50*160*501 - sum(b100_total[,4]))

plot(ct, type1, type="l", ylab="error rate", xlab="minimum minor allele frequency to report", lwd=2.5)
lines(ct, type2, col="red", lwd=2.5)

polygon(c(ct,ct[11:1]), c(type1-ci1,(type1+ci1)[11:1]),
        col=adjustcolor("black",alpha.f=0.15), lty=0)
polygon(c(ct,ct[11:1]), c(type2-ci2,(type2+ci2)[11:1]),
        col=adjustcolor("red",alpha.f=0.15), lty=0)

legend("topright", c("Type I error", "Type II error"), lty=c(1,1), col=c(1,2), lwd=c(2.5,2.5))

dev.off()


