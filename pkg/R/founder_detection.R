#' Detection of founder segments
#'
#' Function to add off-variant-blocks
#' @param dhm haploid SNP-dataset (founders need to be in the first few columns)
#' @param dhm_founder haploid SNP-dataset of the founders
#' @param dhm_founder haploid SNP-dataset of the offspring
#' @param founding Number of founding haplotypes in the population
#' @param plot_overview Generate a plot of the lines chosen in lines_to_plot
#' @param plot_perline Generate an per line plot of estimated founders
#' @param min_condorance Minimum concordance of a segment to be assigned to a founder (default: 0; everything will be assigned)
#' @param multicolor Set to TRUE to indicate regions with multiple potentialy founder haplotypes by use of multiple colors
#' @param err Hidden Markov parameter (chance of error within a haplotype block; default: 0.025
#' @param err2 Hidden Markov parameter (chance of error outside of haplotype blocks; default: 0.1)
#' @param rec Hidden Markov parameter (chance of a recombination event between two markers; default: 0.001)
#' @param ped Pedigree - this will lead to lines only having founders from the lines given ancestors
#' @param weighting_length Weighting factor for length to determine major block (default: 1)
#' @param weighting_size Weighting factor for number of haplotypes in block to determine major block (default: 1)
#' @param min_majorblock minimum of cells in the dataset a block is the biggest covering (default: 5000)
#' @param target_coverage Target Coverage in the blocklist
#' @param big_output If TRUE return not only blocklist but also blockinfo, dataset ect.
#' @return Estimated Founder for each cell in the haploid dataset (dhm)
#' @export
#'
#'

founder_detection <- function(dhm=NULL, dhm_founder = NULL, dhm_off = NULL, founding = NULL, plot_perline=FALSE, plot_overview=TRUE,
                              lines_to_plot=1:50, min_concordance = 0, multicolor=FALSE,
                              err = 0.025, err2 = 0.1, rec = 0.0001, ped = NULL,
                              weighting_length=2,
                              weighting_size=1,
                              target_coverage=NULL,
                              min_majorblock=NULL,
                              big_output = FALSE){

{
  if(length(target_coverage)==0 & length(min_majorblock)==0){
    min_majorblock = 1
  }
  if(length(target_coverage)==1 & length(min_majorblock)==0){
    min_majorblock = 5000
  }

  if(length(dhm)==0){
    dhm = cbind(dhm_founder, dhm_off)
    founding = 1:ncol(dhm_founder)
  }
  if(length(founding)==1){
    founding = 1:founding
  } else if(length(founding)==0){
    stop("Please provide the number of founders via the parameter founding")
  }

  nfounder <- length(founding)
  noff <- ncol(dhm) - nfounder

  recombi_list <- list()

  line_cont <- matrix(0, nrow=nfounder, ncol=ncol(dhm)-nfounder)
  genome_cont <- matrix(0, nrow=nfounder, ncol= nrow(dhm))

  blocklist <- block_calculation(dhm, window_size = 5, min_majorblock = min_majorblock,
                                 subgroups = list(founding, (1:ncol(dhm))[-founding]),
                                 min_per_subgroup = 1,
                                 weighting_length = weighting_length)

  seq_table <- matrix(0, nrow=nrow(dhm), ncol=ncol(dhm)-nfounder)

  t <- coverage_test(blocklist)
  print(mean(t))
  se <- blocklist_startend(blocklist)
  print(mean(se[,2]-se[,1]))



  overl <- numeric(length(blocklist))
  for(index in 1:length(blocklist)){
    overl[index] <- length(intersect(blocklist[[index]][[6]], founding))
  }

  k <- 1
  if(plot_overview){

    X11()
    plot(0,0, ylim=c(1, length(lines_to_plot)+1), xlim=c(1, nrow(dhm)+1), ylab="line", xlab="SNP")
  }

  for(line in (nfounder+1):ncol(dhm)){
    print(line)

    subs <- NULL # Start/end points of haplotype blocks
    foundern <- list() # Potential founders (either of the two haplotypes of a founder)
    for(index in 1:length(blocklist)){
      if(sum(blocklist[[index]][[6]]==line)>0){
        subs <- rbind(subs, c(blocklist[[index]][[2]]$snp, blocklist[[index]][[3]]$snp ))
        foundern[[length(foundern)+1]] <- intersect(blocklist[[index]][[6]], founding)
      }
    }



    # version 2: extend haplotype blocks in case genotypes exactly match (still conservative!)
    if(plot_perline){
      X11()
      plot(0,0, ylim=c(0,8), xlim=c(0,nrow(dhm)), ylab="founder", xlab="SNP")
    }


    avail <- matrix(FALSE, nrow=nrow(dhm), ncol=length(founding))
    for(index in 1:nrow(subs)){
      for(index2 in 1:length(founding)){
        if(sum(foundern[[index]]==index2)>0){
          if(plot_perline){
            polygon(c(subs[index,1], subs[index,2], subs[index,2], subs[index,1]), c(index2-1,index2-1,index2,index2), col=index2, lty=0)
          }


          line_prior <- dhm[subs[index,1]:max(subs[index,1]-50,1),line]
          founder_prior1 <- dhm[subs[index,1]:max(subs[index,1]-50,1),founding[index2]]

          line_down <- dhm[subs[index,2]:min(subs[index,2]+50,nrow(dhm)),line]
          founder_down1 <- dhm[subs[index,2]:min(subs[index,2]+50,nrow(dhm)),founding[index2]]

          con1 <- c(line_prior==founder_prior1 , FALSE)
          firstf1 <- which(con1==FALSE)[1]
          con2 <- c(line_down==founder_down1, FALSE)
          firstf2 <- which(con2==FALSE)[1]

          if(plot_perline){
            polygon(c(subs[index,1]-firstf1+1, subs[index,2]+firstf2-1, subs[index,2]+firstf2-1, subs[index,1]-firstf1+1), c(index2-1,index2-1,index2,index2), col=index2, lty=0)

          }

          avail[(subs[index,1]-firstf1+1):(subs[index,2]+firstf2-2), index2] <- TRUE
        }
      }
    }

    avail <- t(avail)


    if(length(ped)>0){
      avail_parent <- ped[which(ped[,1]==line),-1]
    } else{
      avail_parent <- NULL
    }

    ## Forward - Algorithm:

    forward <- matrix(0, nrow=length(founding), ncol=nrow(dhm))


    forward[avail[,1],1] <- 1-err
    forward[!avail[,1],1] <- err

    if(length(avail_parent)>0){
      forward[-avail_parent,1] <- 0
    }

    forward[,1] <- forward[,1] / sum(forward[,1])

    rec_matrix <- matrix(rec, ncol=nfounder, nrow=nfounder)
    diag(rec_matrix) <- 1 - rec



    mult <- numeric(nfounder)
    for(index in 2:nrow(dhm)){

      mult[avail[,index]] <- 1 - err
      mult[!avail[,index]] <- err

      check1 <- (dhm[index,line] == dhm[index,founding]) & (mult == err)
      mult[check1] <- 1 - err2

      if(length(avail_parent)>0){
        mult[-avail_parent] <- 0
      }

      forward[,index] <- mult * colSums(rec_matrix * forward[,index-1])
      forward[,index] <- forward[,index] / sum(forward[,index])
    }



    ## Backward - Algorithm

    backward <- matrix(0, nrow=length(founding), ncol = nrow(dhm))
    backward[,nrow(dhm)] <- 1 / nfounder
    for(index in (nrow(dhm)-1):1){

      mult[avail[,index+1]] <- 1 - err
      mult[!avail[,index+1]] <- err

      check1 <- dhm[index+1,line] == dhm[index+1,founding]
      mult[check1] <- 1 - err2

      if(length(avail_parent)>0){
        mult[-avail_parent] <- 0
      }

      backward[,index] <- mult * colSums(rec_matrix * backward[,index+1])
      backward[,index] <- backward[,index] / sum(backward[,index])
    }




    ## Forward-Backward

    if(FALSE && rb>0){
      forward[forward<rb] <- rb
      backward[backward<rb] <- rb
    }

    fb <- forward * backward


    stand <- colSums(fb)
    fb <- t(t(fb)/stand)
    seq_line <- numeric(ncol(fb))
    #imp_line <- rep("./.", nrow(map_ref))
    if(plot_perline){
      X11()
      plot(fb[1,], type="l", col=1, lwd=2, ylim=c(0,1))
      for(index in 2:8){
        lines(fb[index,], col=index, lwd=2)
      }
    }

    if(TRUE){
      if(plot_overview) abline(h=k)
      for(index in 1:nrow(dhm)){
        temp <- which(fb[,index] == (max(fb[,index])))
        if(length(temp)>1){
          seq_line[index] <- sample(temp,1)
        } else{
          seq_line[index] <- temp
        }

        if(index>1 && fb[seq_line[index],index]<(fb[seq_line[index-1],index]+0.0001 )){
          seq_line[index] <- seq_line[index-1]
        }
      }
    }

    cont <- matrix(0, nrow=length(founding), ncol = nrow(dhm))
    for(index in 1:ncol(fb)){
      cont[seq_line[index],index] <- 1
    }


    if(min_concordance>0){
      breaks <- unique(c(0,which(diff(seq_line)!=0), length(seq_line)))
      for(index in 1:(length(breaks)-1)){

        cs <- (colSums(dhm[(breaks[index]+1):breaks[index+1], founding, drop=FALSE] == dhm[(breaks[index]+1):breaks[index+1], line]))
        activ <- (which(rowSums(cont[,(breaks[index]+1):breaks[index+1], drop=FALSE])>0))
        ll <- (breaks[index+1]-(breaks[index]))

        if(max(cs[activ]/ll)<min_concordance){
          seq_line[(breaks[index]+1):breaks[index+1]] <- 0
          cont[,(breaks[index]+1):breaks[index+1]] <- 0
        }
      }
    }

    if(multicolor){
      breaks <- unique(c(0,which(diff(seq_line)!=0), length(seq_line)))
      for(index in 1:(length(breaks)-1)){

        cs <- (colSums(dhm[(breaks[index]+1):breaks[index+1], founding, drop=FALSE] == dhm[(breaks[index]+1):breaks[index+1], line]))
        activ <- (which(rowSums(cont[,(breaks[index]+1):breaks[index+1], drop=FALSE])>0))
        ll <- (breaks[index+1]-(breaks[index]))

        if(length(activ)>0){
          if(min(cs[activ]) <= max(cs[-activ])){
            if(length(avail_parent)>0){
              cont[which(cs>=min(cs[activ])),(breaks[index]+1):breaks[index+1]] <- 1 / length(which(cs>=min(cs[activ])))

            } else{
              cont[intersect(avail_parent,which(cs>=min(cs[activ]))),(breaks[index]+1):breaks[index+1]] <- 1 / length(intersect(avail_parent,which(cs>=min(cs[activ]))))
            }

          }
        } else if(sum((cs/ll) > min_concordance)>0){
          if(length(avail_parent)>0){
            cont[which((cs/ll) > min_concordance),(breaks[index]+1):breaks[index+1]] <- 1 / length(which((cs/ll) > min_concordance))
          } else{
            cont[intersect(avail_parent, which((cs/ll) > min_concordance)),(breaks[index]+1):breaks[index+1]] <- 1 / length(intersect(avail_parent, which((cs/ll) > min_concordance)))
          }

        }

      }
    }

    for(index in 1:nrow(dhm)){
      if(plot_overview && !multicolor ) polygon(c(index, index+1, index+1,index), c(k, k, k+1, k+1)-1, col = seq_line[index], lty=0)

      if(plot_overview && multicolor){
        cols <- which(cont[,index]>0)
        ncc <- length(cols)
        if(length(cols)>0){
          for(nc in 1:ncc){
            polygon(c(index, index+1, index+1,index), c(k, k, k+1/ncc, k+1/ncc)-1 + 1/ncc*(nc-1), col = cols[nc], lty=0)
          }
        }


      }
    }

    recombi_event <- c(0,which(diff(seq_line)!=0))



    recombi_list[[k]] <- list()
    recombi_list[[k]][[1]] <- cbind(recombi_event+1, c(recombi_event[-1], nrow(dhm)))
    recombi_list[[k]][[2]] <- list()
    for(index in 1:length(recombi_event)){
      recombi_list[[k]][[2]][[index]] <- which(cont[,recombi_event[index]+1]!=0)
    }

    if(plot_overview) abline(h=k)
    k <- k+1


    seq_table[,k-1] <- seq_line
    line_cont[,line-nfounder] <- rowSums(cont)
    genome_cont <- genome_cont + cont



  }

  if(big_output){
    list(seq_table, line_cont, genome_cont)
  } else{
    return(seq_table)
  }

}
}
