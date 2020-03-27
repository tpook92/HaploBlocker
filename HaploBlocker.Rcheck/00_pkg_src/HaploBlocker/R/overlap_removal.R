#' Function to remove overlapping block segments
#'
#' Function to remove overlapping block segments
#' @param blocklist blocklist
#' @param data window cluster (third output in block_calculation when big_output is set to TRUE)
#' @export
#'

overlap_removal <- function(blocklist=NULL, data=NULL){

  ## Remove all some SNP-extensions
  for(index in 1:length(blocklist)){
    old_start <- blocklist[[index]][[2]]$snp
    old_end <- blocklist[[index]][[3]]$snp

    suppressWarnings(keep <- which(!is.na(as.numeric(blocklist[[index]][[1]]))))
    blocklist[[index]][[1]] <- as.numeric(blocklist[[index]][[1]][keep])

    blocklist[[index]][[2]] <- data[[1]][[min(blocklist[[index]][[1]])]]$start
    blocklist[[index]][[3]] <- data[[1]][[max(blocklist[[index]][[1]])]]$end

    keep2 <- which((old_start:old_end) <= blocklist[[index]][[3]]$snp & (old_start:old_end) >= blocklist[[index]][[2]]$snp)
    blocklist[[index]][[7]]$snp <- blocklist[[index]][[7]]$snp[keep2]
    blocklist[[index]][[7]]$freq <- blocklist[[index]][[7]]$freq[keep2]
  }

  ## Reduce overlap
  t <- coverage_test(blocklist, type="window", max=1000)
  se <- blocklist_startend(blocklist, type="window")
  run <- 0
  while(max(t)>1){
    run <- run +1

    # which position contains the most different blocks
    max_p <- which.max(t)
    indi <- ceiling(max_p / nrow(t))
    window <- max_p - (indi-1)*nrow(t)
    t[window,indi]

    # which blocks are there
    candidates <- which(se[,1] <= window & se[,2]>=window)
    for(index in length(candidates):1){
      if(length(intersect(blocklist[[candidates[index]]][[6]], indi))==0){
        candidates <- candidates[-index]
      }
    }

    new_blocks <- NULL

    main <- blocklist[[candidates[[1]]]][[1]]
    for(index in 2:length(candidates)){
      main <- base::intersect(main, blocklist[[candidates[[index]]]][[1]])
    }
    main <- suppressWarnings(as.numeric(main))
    main <- main[!is.na(main)]

    nodes <- main
    if(length(main)>0){
      start <- data[[1]][[min(main)]][[1]]
      end <- data[[1]][[max(main)]][[2]]
      infostuff <- NULL
      for(index2 in nodes){
        infostuff <- c(infostuff,    data[[1]][[index2]][[4]] )
      }
      indis <- NULL
      for(index in candidates){
        indis <- c(indis, blocklist[[index]][[6]])
      }
      indis <- sort(unique(indis))


      new_blocks[[length(new_blocks)+1]] <- list(nodes,start,end,infostuff,length(indis), indis)
    } else{
      seg <- mean(max(se[candidates,1]), min(se[candidates,2]))
      main <- suppressWarnings(as.numeric(blocklist[[candidates[[1]]]][[1]]))
      main <- main[!is.na(main)]
      pos <- numeric(length(main))
      for(index in 1:length(main)){
        pos[index] <- data[[1]][[main[index]]]$start$window
      }
      main <- suppressWarnings(c(max(main[pos<seg]), min(main[pos>seg])))
      main <- main[main>0]
      main <- main[main<Inf]
      main <- mean(main, na.rm =TRUE)
    }


    upstream <- list()
    downstream <- list()

    for(index in 1:length(candidates)){
      upstream[[index]] <- sort(as.numeric(blocklist[[candidates[[index]]]][[1]][blocklist[[candidates[[index]]]][[1]]<min(main)]), decreasing = TRUE)
      downstream[[index]] <- blocklist[[candidates[[index]]]][[1]][blocklist[[candidates[[index]]]][[1]]>max(main)]
      downstream[[index]] <- suppressWarnings(as.numeric(downstream[[index]]))
      downstream[[index]] <- downstream[[index]][!is.na(downstream[[index]])]
    }

    #upstream blocks
    lblock <- llength(upstream)
    while(max(lblock)>0){
      active_block <- which(lblock>0)
      seg_length <- min(lblock[active_block])
      variants <- matrix(0, nrow=length(active_block), ncol=seg_length)
      active_indi <- NULL
      for(index in 1:length(active_block)){
        variants[index,] <- upstream[[active_block[index]]][1:seg_length]
        active_indi <- c(active_indi, blocklist[[candidates[[active_block[index]]]]][[6]])
      }
      active_indi <- sort(unique(active_indi))
      block_variants <- unique(variants)

      for(index in 1:nrow(block_variants)){

        nodes <- sort(block_variants[index,])
        start <- data[[1]][[min(block_variants)]][[1]]
        end <- data[[1]][[max(block_variants)]][[2]]
        infostuff <- NULL
        indis <- active_indi
        for(index2 in nodes){
          infostuff <- c(infostuff,    data[[1]][[index2]][[4]] )
          indis <- intersect(indis, data[[1]][[index2]][[5]] )
        }


        new_blocks[[length(new_blocks)+1]] <- list(nodes,start,end,infostuff,length(indis), indis)
      }
      for(index in active_block){
        upstream[[index]] <- upstream[[index]][-(1:seg_length)]
      }
      lblock <- llength(upstream)

    }

    # downstream blocks
    lblock <- llength(downstream)
    while(max(lblock)>0){
      active_block <- which(lblock>0)
      seg_length <- min(lblock[active_block])
      variants <- matrix(0, nrow=length(active_block), ncol=seg_length)
      active_indi <- NULL
      for(index in 1:length(active_block)){
        variants[index,] <- downstream[[active_block[index]]][1:seg_length]
        active_indi <- c(active_indi, blocklist[[candidates[[active_block[index]]]]][[6]])
      }
      active_indi <- sort(unique(active_indi))
      block_variants <- unique(variants)

      for(index in 1:nrow(block_variants)){

        nodes <- sort(block_variants[index,])
        start <- data[[1]][[min(block_variants)]][[1]]
        end <- data[[1]][[max(block_variants)]][[2]]
        infostuff <- NULL
        indis <- active_indi
        for(index2 in nodes){
          infostuff <- c(infostuff,    data[[1]][[index2]][[4]] )
          indis <- intersect(indis, data[[1]][[index2]][[5]] )
        }


        new_blocks[[length(new_blocks)+1]] <- list(nodes,start,end,infostuff,length(indis), indis)
      }
      for(index in active_block){
        downstream[[index]] <- downstream[[index]][-(1:seg_length)]
      }
      lblock <- llength(downstream)

    }

    for(index in sort(candidates, decreasing = TRUE)){
      blocklist[[index]] <- NULL
    }
    blocklist <- c(blocklist, new_blocks)

    t <- coverage_test(blocklist, type="window", max=1000)
    se <- blocklist_startend(blocklist, type="window")

  }

  # SNP extension without Overlap

  new_order <- sort(se[,1], index.return=TRUE)$ix
  old_blocklist <- blocklist
  blocklist <- NULL
  for(index in 1:length(old_blocklist)){
    blocklist[[index]] <- old_blocklist[[new_order[index]]]
  }
  return(blocklist)

}
