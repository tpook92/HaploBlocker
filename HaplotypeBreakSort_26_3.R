###### HaploBlocker - sorting pipeline
##############################################
##  Requirement: HaploBlocker 1.5.4         ##
##  Requirement: RandomFieldsUtils 0.5.9    ##
##############################################

####### Import parameter settings
# Arg components are supposed to be:
# 1. path of the json-file
# 2. number of bins
# 3. Number of Break points to identify
# 4. MAF-filter (If not set apply no filtering)
# 5. Window size in HB (blocks will be larger!)
# 6. No overlap between blocks (sorting itself seems to be better with overlapping blocks)
# 7. Generate HaploBlocker visualization when TRUE
# 8. Indicate regions with present paths and no blocks
# 9. Merge path of a specimens to a joint path for the sorting
# 10. Minimum size of the haplotype blocks (bins covered by each haplotype block)

args <- commandArgs(TRUE)
# Example args for Yeast 1010
# args <- c("C:/Users/pook/Desktop/chrI.1000.w100.wo_links.json", 3, 0.002, 5, T, T, T, F, 500)
# Rscript --vanilla /mnt/vol1/software/HaploBlocker/HaplotypeBreakSort.R /mnt/vol1/tg/40_odgi/10_yeast/30_sebastian/chrI.1000.w100.wo_links.json 3 0.002 5 T T T F 500
args <- c("C:/Users/pook/Desktop/athaliana_10_02/sebastian.Athaliana.all.50000.w1000.json", 10, 0.002, 5, T, T, T, F, 100)
args <- c("C:/Users/pook/Desktop/graph-genome-data/chrk_ath_12samples.w1000.json", 10, 0.002, 5, T, T, T, F, 100)

args <- c("sebastianOLD.Athaliana.all.50000.w1000.json", 10, 0.002, 5, T, T, T, F, 50)

args <- c("C:/Users/pook/Desktop/graph-genome-data/run1.B1phi1.i1.seqwish.w100.json", 10, 0.002, 5, T, T, T, F, 100)

args <- c("C:/Users/pook/Desktop/graph-genome-data/chrk_ath_80samples_10kb.w1000.json", 10, 0.002, 5, F, T, T, F, 100)

args <- c("C:/Users/pook/Desktop/graph-genome-data/chrk_ath_80samples_10kb_Chr1.w100000.json", 10, 0.002, 5, F, T, T, F, 100)

#args <- c("C:/Users/pook/Desktop/run1.B1phi1.i1.seqwish.w100.json", 10, 0.002, 5, T, T, T, F, 500)
jsonpath <- args[1]
nbreak <- as.numeric(args[2]) # number of breaks we want to have
maf <- as.numeric(args[3])
window_size <- as.numeric(args[4]) # HB window size (Smaller than blocks we are looking for!)
overlap_remove <- as.logical(args[5])
generate_plot <- as.logical(args[6]) # generate a png file with the HB graph
show_black <- as.logical(args[7]) # Display present regions in no block in black
merge_path <- as.logical(args[8]) # Merge all path into a joint path
min_majorblock <- as.numeric(args[9]) # Min number of SNPs each block is covering
nbins <- as.numeric(args[10]) # Number of bins // saves me time to compute if provided but optional
## Current problem: When there are multiple path for a specimen target_coverage is not a good method
## because a high share of coverage is not expected
## To-do: paralization of new haploblocker code

# Some parameters initialize with some random values - will work on automated values for those
min_haplotypes_per_block <- NULL # Default based on number of individuals in panel (>500 -> 5, >25 --> 3, else 2)
bandwidth <- 25 # some smoothing for highly linked potential break points
n_colors <- 20 # number of colors to use visualization (install RcolorBrewer)
min_ori <- 5 # Minimum number of blocks to consider in the ordering process
min_distance <- NULL # Use this to hard code minimum distance between break points

library(HaploBlocker) # Obvious
library(jsonlite) # Write the output json
library(RColorBrewer) # Color scheme for HB visualization

if(length(nbins)==0 || is.na(nbins)){
  path_set1 <- scan(jsonpath, what="character", nlines = 1)
  nbins <- as.numeric(strsplit(path_set1[4], split="}")) /   as.numeric(strsplit(path_set1[3], split=",")[[1]][1])
  nbins <- ceiling(nbins) + 1
}

if(length(nbins)==0 || is.na(nbins)){
  path_set1 <- scan(jsonpath, what="character")
  npath <- 0
  while(substr(path_set1[[length(path_set1)-npath]], start=1, stop=5)=='{\"pat'){
    npath <- npath + 1
  }
  nbins <- length(path_set1)-npath
}

path_set <- scan(jsonpath, what="character", skip=nbins)
indi <- length(path_set)
exc <- numeric(indi)
chr <- numeric(indi)

# Reading out the paths - this is a slow R implementation using strsplit.
# Can be improved when it becomes a computationally relevant part of the script
paths <- list()
coverages <- list()
for(index in 1:length(path_set)){
  print(index)
  path1 <- strsplit(path_set[[index]], split="[", fixed=TRUE)
  if(length(path1)>nbins){
    a <- 0
    for(index in 3:nbins){
      if(length(path1[[index]])==2){
        a <- index -1
        break
      }
    }
    path1[[1]] <- path1[[1]][1:a]
  }
  path1 <- strsplit(path1[[1]], split=",")
  vec <- coverage <- numeric(min(length(path1)-3, nbins))

  for(index2 in 3:(length(path1)-1)){
    if(length(path1[[index2]])>=4 && path1[[index2]][1]!=0){
      vec[index2-2] <- path1[[index2]][1]
      coverage[index2-2] <- path1[[index2]][2]
    } else{
      break
    }
  }

  temp1 <- strsplit(path1[[1]], split='\"', fixed=TRUE)
  temp1 <- strsplit(temp1[[1]], split="_Chr", fixed=TRUE)
  exc[index] <- temp1[[length(temp1)]][1]
  chr[index] <- temp1[[length(temp1)]][2]
  if(is.na(chr[index])){
    temp1 <- strsplit(exc[index], split=".Chr", fixed=TRUE)
    exc[index] <- temp1[[1]][1]
    chr[index] <- temp1[[1]][2]
  }
  paths[[index]] <- as.numeric(vec)
  coverages[[index]] <- as.numeric(coverage)
}

nsnp <- max(unlist(paths))

# Generate input dataset for haploblocker
dataset1 <- matrix(-9, nrow=nsnp, ncol=indi)
for(index in 1:indi){
  dataset1[paths[[index]], index] <- coverages[[index]][paths[[index]]!=0]
}
markers_kept <- which(rowSums(dataset1!=(-9))>(indi*maf))
dataset2 <- dataset1[markers_kept,]



## Dataset coding
if(TRUE){
  dataset <- dataset2
  dataset[dataset!=(-9)] <- 1
} else if(TRUE){
  dataset <- dataset2
  dataset[dataset2>=0.1] <- 1
  dataset[dataset2<0.1] <- (-9)
} else if(TRUE){
  dataset <- dataset2
  dataset[dataset2< 0.1 & dataset2 >0] <- 1
  dataset[dataset2< 0.4 & dataset2 >=0.1] <- 2
  dataset[dataset2< 0.8 & dataset2 >=0.4] <- 3
  dataset[dataset2>= 0.8] <- 4
} else if(TRUE){
  dataset3 <- dataset2
  dataset3[dataset3<0] <- 0
  dataset3[dataset3>2] <- 2
  p <- rowMeans(dataset3)
  dataset <- dataset3
  dataset[dataset3<(0.8*p)] <- 1
  dataset[dataset3<(1.2*p) & dataset3>=(0.8*p)] <- 2
  dataset[dataset3>=(1.2*p)] <- 3
  dataset[dataset3<=0.1] <- -9
}

if(length(min_haplotypes_per_block)==0){
  if(indi>500){
    min_haplotypes_per_block <- 5
  } else if(indi>100){
    min_haplotypes_per_block <- 3
  } else{
    min_haplotypes_per_block <- 2
  }
}

if(merge_path){
  exc1 <- unique(exc)

  dataset2 <- matrix(0L, nrow=nsnp, ncol=length(exc1))

  for(index in 1:length(exc1)){
    dataset2[,index] <- rowSums(dataset1[,which(exc==exc1[index]), drop=FALSE])
  }
  dataset <- dataset2[rowSums(dataset1)>(indi*maf),]
}

## Apply HaploBlocker
system.time({
  #dataset <- dataset[rowMeans(dataset)>(-8.9),]
  temp11 <- block_calculation(dataset[1:25000,], node_min=min_haplotypes_per_block,
                              window_size = window_size,
                              merging_error = if(window_size==5){1}else{2},
                              weighting_length = 2,
                              weighting_size = 1,
                              edge_min = min_haplotypes_per_block,
                              min_majorblock = min_majorblock,
                              deletion_count=TRUE,
                              big_output = TRUE,
                              overlap_remove = overlap_remove,
                              off_lines = 0, parallel_window = 5000,
                              window_overlap = 0,
                              window_cores=4, gap=0)
})
temp11 <- block_calculation(dataset, node_min=min_haplotypes_per_block,
                            window_size = window_size,
                            merging_error = if(window_size==5){1}else{2},
                            weighting_length = 2,
                            weighting_size = 1,
                            edge_min = min_haplotypes_per_block,
                            min_majorblock = min_majorblock,
                            deletion_count=TRUE,
                            big_output = TRUE,
                            overlap_remove = overlap_remove,
                            off_lines = 0, gap=0)

'#

temp11 <- block_calculation(dataset, node_min=min_haplotypes_per_block,
window_size = window_size,
edge_min = min_haplotypes_per_block,
min_majorblock = min_majorblock,
deletion_count=TRUE,
big_output = TRUE,
overlap_remove = overlap_remove,
parallel_window=10000,
window_overlap=1000,
window_cores=4)

'#



blocklist <- temp11[[1]]
present_data <- temp11[[8]]

for(index in 1:length(blocklist)){
  keep <- which(blocklist[[index]][[7]]$snp!=(-9))
  if(length(keep)>0){
    kstart <- min(keep)
    kend <- max(keep)
    blocklist[[index]][[2]]$snp <- blocklist[[index]][[2]]$snp + kstart - 1
    blocklist[[index]][[3]]$snp <- kend - kstart + blocklist[[index]][[2]]$snp
    blocklist[[index]][[2]]$window <- ceiling(blocklist[[index]][[2]]$snp / window_size)
    blocklist[[index]][[3]]$window <- ceiling(blocklist[[index]][[3]]$snp / window_size)
    blocklist[[index]][[7]]$snp <- blocklist[[index]][[7]]$snp[kstart:kend]
    blocklist[[index]][[7]]$freq <- blocklist[[index]][[7]]$freq[kstart:kend]
  }
}

## Some filtering to remove blocks with no path + make sure start/end point are present path
for(index in length(blocklist):1){
  if((sum(blocklist[[index]][[7]]$snp!=(-9))/ length(blocklist[[index]][[7]]$snp))<0.03){
    blocklist[[index]] <- NULL
  }
}

if(TRUE){
  t <- t(coverage_test(blocklist, type="window"))
  t <- cbind(t, matrix(0, ncol=ncol(present_data) -ncol(t), nrow=nrow(t) ))
  cat(paste0( round(sum(present_data[,1:ncol(t)] * t) / sum(present_data)*100, digits=2), "% of the existing variation is in haplotype blocks!\n"))
}

# Just some evaluations to make sure results are resonable
if(FALSE){

  t <- coverage_test(blocklist, type="window")
  t[1:72,51]==1
  present_data[51,1:72]
  se <- blocklist_startend(blocklist)
  mean(se[,2]-se[,1])
  median(se[,2]-se[,1])
  share_missing <- numeric(length(blocklist))
  for(index in 1:length(blocklist)){
    share_missing[index] <- sum(blocklist[[index]][[7]]$snp!=(-9)) / length(blocklist[[index]][[7]]$snp)
  }
  mean(coverage_test(blocklist, type="snp"))
  mean(share_missing)
}


## How many block have break points in a location
## For how specimens have blocks the span over the break point
se <- blocklist_startend(blocklist)
end_block  <- sort(unique(c(0,se[,1]-1, se[,2])))[-1]
freq_end <- freq_end_jump <- numeric(length(end_block))

{
  vars <- list()
  for(index in 1:length(blocklist)){
    count <- which((blocklist[[index]][[2]]$snp-1) == end_block)
    if(length(count)>0){
      if(length(vars)>=count){
        vars[[count]] <- unique(c(vars[[count]], blocklist[[index]][[6]]))
      } else{
        vars[[count]] <- blocklist[[index]][[6]]
      }
    }

    count <- which((blocklist[[index]][[3]]$snp) == end_block)
    if(length(count)>0){
      if(length(vars)>=count){
        vars[[count]] <- unique(c(vars[[count]], blocklist[[index]][[6]]))
      } else{
        vars[[count]] <- blocklist[[index]][[6]]
      }
    }


  }
  for(index in 1:length(vars)){
    freq_end[index] <- length(vars[[index]])
  }

  # How many individuals are in the same block before each block-end?
  for(index in 1:length(end_block)){
    set <- NULL
    candidate_block <- which(((se[,1]+window_size)<=end_block[index]) & ((se[,2]-window_size)>(end_block[index]+1)))
    for(index2 in candidate_block){
      set <- unique(c(set, blocklist[[index2]][[6]]))
    }
    freq_end_jump[index] <- length(set)
  }
}

#  plot(end_block, freq_end)
#  plot(end_block, freq_end_jump)
freq1 <- freq2 <- numeric(max(se))
for(index in 1:length(end_block)){
  freq1[end_block[index]] <- freq_end[index]
  freq2[end_block[index]] <- freq_end_jump[index]
}
#  plot(freq1)
#  plot(freq1 * (1-freq2/501))

###########################################
#### All counting is done #################
###########################################
### Start the visualization ###############
###########################################

nsize <- blocklist_size(blocklist)

## some plot_block defaults
intensity <- 0.7
min_to_plot <- 2
perform_order <-TRUE
add_sort <- TRUE
xlim <- NULL
max_step <- 500

if(length(min_distance)==0){
  min_distance <- ceiling(max(se) / nbreak / 2 / 4)
}
### Derive Break points
{
  breaks <- NULL
  if(max(freq2/indi)<0.5){
    freq22 <- freq2*indi/max(c(freq2,1))
  } else{
    freq22 <- freq2
  }

  r1 <- freq1 * (1-freq22/indi)
  bandwidth <- min(bandwidth, min_distance)
  a <- ksmooth(1:max(se), r1, bandwidth = bandwidth, kernel="normal")

  #plot(a, type = "l")
  if(min_distance<Inf){
    a$y[1:min_distance] <- 0 # no break points at start/end
    a$y[(length(a$y)-min_distance) : length(a$y)] <- 0
  }

  freq11 <- freq1
  if(nbreak>0){
    for(index in 1:nbreak){
      candidate <- which.max(a$y)
      if(candidate==1 && a$y[1]==0){
        break()
      }
      new_break <- which.max(freq11[(candidate - bandwidth):(candidate + bandwidth)]) + candidate - bandwidth - 1
      nb2 <- new_break
      new_break <- end_block[which.min(abs(new_break - end_block))]
      breaks <- c(breaks, new_break)
      a$y[(min(new_break, nb2)-min_distance) : (max(new_break, nb2) + min_distance)] <- 0
      freq11[(new_break-min_distance) : (new_break + min_distance)] <- 0
    }
  }

  breaks <- sort(breaks)
}


### Derive ordering of haplotypes per window
{
  windows <- c(1, breaks, max(se))

  present_score <- matrix(0, nrow=indi, ncol=length(windows)-1)
  for(win in 1:(length(windows)-1)){
    present_score[,win] <- rowMeans(present_data[,floor(windows[win]/window_size):floor((windows[win+1])/window_size), drop=FALSE])
  }
  anchors <- numeric(length(windows)-1)
  order_window <- list()
  bl_size <- blocklist_size(blocklist)
  t_cov1 <- t(coverage_test(blocklist, type="snp"))
  t_cov <- numeric(ncol(t_cov1))
  for(win in 1:(length(windows)-1)){
    order <- NULL
    snp_ori <- ceiling(mean(windows[win:(win+1)]))

    t_cov[windows[win]:(windows[win+1]-1)] <- colMeans(t_cov1[,windows[win]:(windows[win+1]-1)] * present_score[,win]) ### SCORE CALCULATION!
    if(max(t_cov[windows[win]:(windows[win+1]-1)]) > (t_cov[snp_ori]*1.2)){
      snp_ori <- which.max(t_cov[windows[win]:(windows[win+1]-1)]) + windows[win] -1
      if(snp_ori==windows[win]){
        snp_ori <- snp_ori +1
      }
      if(snp_ori==windows[win+1]){
        snp_ori <- snp_ori -1
      }
    }
    anchors[win] <- snp_ori
    orientation <- which(((se[,1]<= snp_ori)+ (se[,2]>= snp_ori))==2)
    stretch <- 0
    while(length(orientation)<min_ori){
      stretch <- stretch + 10
      orientation <- which(((se[,1]<= (snp_ori+stretch))+ (se[,2]>= (snp_ori-stretch)))==2)
    }
    non_included <- se[orientation,1]>windows[win+1] | se[orientation,2]<windows[win]
    orientation <- orientation[!non_included]
    take <- which.max(nsize[orientation])
    sorted <- orientation[take]
    orientation <- orientation[-take]
    while(length(orientation)>0){
      overlap <- numeric(length(orientation))
      actives <- blocklist[[sorted[length(sorted)]]][[6]]
      i1 <- 1
      for(index in orientation){
        nextone <- blocklist[[index]][[6]]
        overlap[i1] <- length(intersect(nextone,actives))/ length(nextone)
        i1 <- i1 + 1
      }
      sorted <- c(sorted, orientation[which.max(overlap)[1]])
      orientation <- orientation[-which.max(overlap)[1]]
    }
    orientation <- sorted

    if(TRUE){
      orientation2 <- which(((se[,1]>= windows[win])+ (se[,2]<= windows[win+1]))==2)
      orientation2 <- orientation2[sort(bl_size[orientation2], index.return=TRUE, decreasing=TRUE)$ix]
      orientation <- c(orientation, orientation2)
      orientation <- orientation[!duplicated(orientation)]
    }
    if(TRUE){
      orientation3 <- which((se[,1]>=windows[win] & se[,1]<windows[win+1]) | (se[,2]>=windows[win] & se[,2]<windows[win+1]))
      orientation3 <- orientation3[sort(bl_size[orientation3], index.return=TRUE, decreasing=TRUE)$ix]
      orientation <- c(orientation, orientation3)
      orientation <- orientation[!duplicated(orientation)]
    }


    order <- NULL
    orientation <- unique(c(orientation, 0))
    if(perform_order){
      for(index in orientation){
        pl <- length(order)
        if(index==0){

          if(TRUE){
            orientation4 <- sort(present_score[,win], index.return=TRUE, decreasing=TRUE)$ix
            remainder <-orientation4
          } else{
            remainder <- 1:indi
          }


          order <- unique(c(order,remainder))
          index <- ceiling(median(orientation[-length(orientation)]))
        } else{
          order <- unique(c(order,blocklist[[index]][[6]]))
        }

        if(add_sort==TRUE && length(order)>pl){
          if(pl>0){
            added <- order[-(1:pl)]
          } else{
            added <- order
          }
          groups <- list()
          groups[[1]] <- added
          for(step in 1:max_step){

            if(length(groups)<length(added)){
              if((index-step) >0 && blocklist[[index-step]][[3]]$snp > windows[win] && length(intersect(blocklist[[index-step]][[6]], added))>0){
                new_groups <- list()
                for(group in 1:length(groups)){
                  same <- base::intersect(groups[[group]], blocklist[[index-step]][[6]])
                  rest <- base::intersect(groups[[group]], (1:indi)[-blocklist[[index-step]][[6]]])
                  if(length(same)>0){
                    new_groups[[length(new_groups)+1]] <- same
                  }
                  if(length(rest)>0){
                    new_groups[[length(new_groups)+1]] <- rest
                  }
                }
                groups <- new_groups
              }
              if((index+step) <=length(blocklist) && blocklist[[index+step]][[2]]$snp < windows[win+1] && length(intersect(blocklist[[index+step]][[6]], added))>0){
                new_groups <- list()
                for(group in 1:length(groups)){
                  same <- base::intersect(groups[[group]], blocklist[[index+step]][[6]])
                  rest <- base::intersect(groups[[group]], (1:indi)[-blocklist[[index+step]][[6]]])
                  if(length(same)>0){
                    new_groups[[length(new_groups)+1]] <- same
                  }
                  if(length(rest)>0){
                    new_groups[[length(new_groups)+1]] <- rest
                  }
                }
                groups <- new_groups
              }
            }
          }
          if(pl>0){
            order <- c(order[(1:pl)],unlist(groups))
          } else{
            order <- c(unlist(groups))
          }

        }

      }

      order_window[[win]] <- order
    } else{
      order <-  1:indi
      order_window[[win]] <- order
    }

  }
}

{
  if(show_black){
    back_order <- list()
    for(win in 1:(length(windows)-1)){
      back <- numeric(length(order_window[[win]]))
      for(index2 in 1:length(back)){
        back[index2] <- which(order_window[[win]]==index2)
      }
      back_order[[win]] <- back
    }
  }
}


### Actual plotting
jsonpath_short <- substr(jsonpath, start=1, stop=(nchar(jsonpath)-5))
if(generate_plot){

  order_window[[1]] <- 1:ncol(dataset)
  back_order[[1]] <- 1:ncol(dataset)
  #png(file=paste0(jsonpath, ".png"), width=2250, height= 960, res=300)
  par(mar=c(3.1,2.1,0.1,0.1))
  if(length(xlim)>0){
    plot(0,-1000,ylim=c(0,length(order)), ylab="haplotype", xlab="SNP",
         cex.axis=1, cex.lab=1, xlim=xlim)
  } else{
    plot(0,-1000,ylim=c(0,length(order)), xlim=c(1,max(se)), ylab="haplotype", xlab="SNP",
         cex.axis=1, cex.lab=1)
  }

  if(requireNamespace("RColorBrewer")){
    qual_col_pals = RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == 'qual',]
    col_vector = unlist(mapply(RColorBrewer::brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    used_color <- sample(col_vector, n_colors)
  } else{
    used_color <- 1:8
    n_colors <- 8
  }
  activ_colors <- numeric(length(blocklist))
  for(win in 1:(length(windows)-1)){
    activ_block <- which((se[,1] < windows[win+1]) & (se[,2] > windows[win]))

    for(index in activ_block){
      overlap <- duplicated(c(blocklist[[index]][[6]], order_window[[win]]))[-(1:blocklist[[index]][[5]])]
      if(sum(overlap) >= min_to_plot){
        if(activ_colors[index]==0){
          taken <- base::intersect(which(se[,1]<=se[index,2]), which(se[,2]>=se[index,1]))
          block_color <- sort(unique(c(0,activ_colors[taken])))
          if(length(block_color)>1 && length(block_color)<=n_colors){
            activ_colors[index] <- sample((1:n_colors)[-block_color],1)
          } else{
            activ_colors[index] <- sample((1:n_colors),1)
          }
        }

        if(TRUE){
          poly <- sort(which(overlap))
          takes <- c(TRUE,diff(poly)>1)
          to_plot <- poly[takes]
          size <- diff(which(takes))
          size <- c(size, length(poly)-sum(size))
          temp1 <- 1
          for(index2 in to_plot){
            polygon(c(max(se[index,1], windows[win]), min(se[index,2], windows[win+1]),
                      min(se[index,2], windows[win+1]), max(se[index,1], windows[win])), index2-c(1,1,1-size[temp1],1-size[temp1]),
                    col=adjustcolor(used_color[activ_colors[index]],alpha.f=intensity), lty=0)
            temp1 <- temp1 +1


          }

        }

      }
    }

    if(show_black){
      start <- ceiling(windows[win]/window_size)
      activ_block <- as.numeric(present_data[, (start):(windows[win+1]/window_size), drop=FALSE]>0) * (t[,(start):(windows[win+1]/window_size), drop=FALSE]==0)
      for(index in 1:indi){
        bbox <- which(activ_block[index,]>0)
        takes <- c(TRUE,diff(bbox)>1)
        to_plot <- bbox[takes]
        size <- diff(which(takes))
        size <- c(size, length(bbox)-sum(size))
        temp1 <- 1
        for(index2 in to_plot){

          polygon(c(index2*window_size-window_size, (index2+size[temp1]-1)*window_size,
                    (index2+size[temp1]-1)*window_size, index2*window_size-window_size)+ (start-1)*window_size, back_order[[win]][index]-c(1,1,0,0),
                  col="black", lty=0)
          temp1 <- temp1 + 1
        }
      }

    }
  }

  abline(v=breaks, lwd=0.4, col="red")
  abline(v=anchors, lwd=0.2, col="black")
  #dev.off()
}

mean(se[,2]-se[,1])
nb <- numeric(max(se))
for(index in 1:nrow(se)){
  nb[se[index,1]:se[index,2]] <- nb[se[index,1]:se[index,2]]+1
}
mean(nb)
plot(nb)

breakpoints <- list()
test <- list()
for(index in 1:(length(windows)-1)){
  breakpoints[[index]] <- list()
  if(index==1){
    breakpoints[[index]]$first_bin <- 1
  } else{
    breakpoints[[index]]$first_bin <- markers_kept[windows[index]]
  }
  if(index==(length(windows)-1)){
    breakpoints[[index]]$last_bin <- nsnp
  } else{
    breakpoints[[index]]$last_bin <- markers_kept[windows[index+1]]-1
  }


  breakpoints[[index]]$row_order <- paste0(exc[order_window[[index]]], ".chr", chr[order_window[[index]]])

}

test$breakpoints <- breakpoints

write_json(test, path=paste0(jsonpath_short, "_breaks.json"))
