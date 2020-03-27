pkgname <- "HaploBlocker"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
base::assign(".ExTimings", "HaploBlocker-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('HaploBlocker')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("coding")
### * coding

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: coding
### Title: Coding
### Aliases: coding codeSNPs decodeSNPs
### Keywords: models

### ** Examples

animals <- 8
snps <- 25
ACGT <- c("A", "C", "G", "T")
fixcoding(c(ACGT, "ANY"))
(M <- matrix(nc=animals, sample(c(ACGT, "-", "+"), animals*snps, replace=TRUE)))
(CM <- codeSNPs(M))
(decCM <- decodeSNPs(CM))
stopifnot(all(M == decCM || ((M == "-" || M == "+") && decCM == "@") ))




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("coding", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("colSumsEqualSNPs")
### * colSumsEqualSNPs

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: colSumsEqualSNPs
### Title: Columwise comparison of two genetic matrices
### Aliases: colSumsEqualSNPs
### Keywords: misc

### ** Examples

fixcoding(0:1)
require(RandomFieldsUtils)
for (i in 1:10) {
  animals <- sample(100, 1)
  snps <- sample(30, 1)
  v_snps <- sample(snps, 1)
  start <- sample(snps - v_snps + 1, 1)

  M <- matrix(nc=animals, sample(c(0,1), animals * snps, replace = TRUE))
  V <- sample(c(0,1), v_snps, replace = TRUE)
  Vext <- c(rep(NA, start-1), V, rep(NA, snps - v_snps - start + 1))
  stopifnot(length(Vext) == nrow(M))
  CM <- codeSNPs(M)
  e <- colSumsEqualSNPs(CM, start, V)
  Print(animals, snps, v_snps, start, e)
  confirm(e, colSums(M == Vext, na.rm=TRUE))
}




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("colSumsEqualSNPs", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("factorSNPs")
### * factorSNPs

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: factorSNPs
### Title: Determines which codes are used in a given part of a SNP
###   sequence
### Aliases: factorSNPs
### Keywords: misc

### ** Examples

require(RandomFieldsUtils)

fixcoding(0:1)
for (i in 1:10) {
  cat(i, "")
  animals <- sample(500, 1)
  snps <- sample(70, 1)
  v_snps <- sample(snps - 1, 1)
  start <- sample(snps - v_snps, 1)
  end <- sample(start:snps, 1)
  Print(animals, snps, start, end)
  
  M <- matrix(nc=animals, sample(c(0,1), animals * snps, replace = TRUE))
  M <- cbind(M, M)
  M <- M[, sample(ncol(M)), drop=FALSE]
  ##  print(M)
  CM <- codeSNPs(M)
  f <- factorSNPs(CM, start, end)
  ff <- unique(M[start:end, ,drop=FALSE], MARGIN = 2)
  ##
  stopifnot(ncol(ff) == length(attr(f, "where.to.find")))
  stopifnot(ncol(ff) == length(attr(f, "where.to.find")))
  stopifnot(all(ff == M[start:end, attr(f, "where.to.find")]))
}



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("factorSNPs", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("fixcoding")
### * fixcoding

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: fixcoding
### Title: Determination of Coding
### Aliases: fixcoding
### Keywords: misc

### ** Examples

fixcoding(0:1)

animals <- 8
snps <- 25
(M <- matrix(nc=animals, sample(0:1, animals*snps, replace=TRUE)))
(CM <- codeSNPs(M))
(decCM <- decodeSNPs(CM))
confirm(M, decCM)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("fixcoding", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("intersect")
### * intersect

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: intersect
### Title: Intersection of sets
### Aliases: intersect
### Keywords: misc

### ** Examples

(a <- as.integer(sort(sample(20, 6))))
(b <- as.integer(sort(sample(20,7))))
confirm(intersect(a,b), base::intersect(a,b))

a <- as.integer(sort(sample(2000, 600)))
b <- as.integer(sort(sample(2000,700)))
confirm(intersect(a,b), base::intersect(a,b))
print(system.time(for (i in 1:10^5) base::intersect(a,b)))
print(system.time(for (i in 1:10^5) intersect(a,b))) ## factor 10 faster




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("intersect", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
