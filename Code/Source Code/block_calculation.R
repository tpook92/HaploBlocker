#' Mainfunction to calculate haplotype blocks
#'
#' Function to generate haplotype blocks from haploid data
#' @param dhm haploid SNP-dataset
#' @param window_size size of each window in the algorithm (default: 20)
#' @param merging_error number of allowed errors per block (default: 1)
#' @param window_sequence sequence of predefined windows (default: NULL ;per row: start$snp, end$snp, length, length - merging_error, start$bp, end$bp)
#' @param bp_map vector of positions for each SNP in bp (default: NULL - all 0)
#' @param window_anchor_gens matrix to constructed window_sequence base on start/end points in bp (e.g. gen regions, per row: start, end)
#' @param at_least_one If TRUE no allowed merging errors in windows of size 1
#' @param minimum.blocksize minimum number of haplotypes per block (default: 5)
#' @param gap remove haplotypes with knoten over minimum.blocksize in less windows (default: 10)
#' @param mindestanteil minimum percentage of transition to the same block for extension (default: 0.95, step III)
#' @param off_lines minimum number of haplotypes to looose in the creation of a bigger block (default: 5)
#' @param mindestaehnlichkeit minimum rate of the same SNPs to be added to the block (default: 0.99)
#' @param merge_closeblock If TRUE merge adjecting block with similar haplotypes (default: FALSE)
#' @param max_diff_l maximum number of windows with different haplotypes inbetween (default: 1)
#' @param max_diff_i maximum number of individuals in only one of the two blocks (default: 1)
#' @param block_min_count minimum of positions in the dataset a block is the biggest covering (default: 50)
#' @param consider_knoten Use knots to identify blocks (default: TRUE)
#' @param consider_trans Use edges between knots to identify blocks (default: TRUE)
#' @param trans_min minimum number of haplotypes per transition to use in consider_trans (default: 5)
#' @param subgroups possible subgroups to consider in the block identification process (default: NULL - list(1:indi))
#' @param min_per_subgroup minimum number of haplotypes per block per subgroup (default: 0)
#' @param consider_all If TRUE always haplotypes which are not in the knot to be in a generated block
#' @param save_allblock If TRUE keep all haplotypes with all windows according to a block (even under mindestaehnlichkeit)
#' @param block_extending If TRUE use the window-extending algorithm (step V)
#' @param max_extending_diff Maximum number of windows with different realisation in the block-extending-algorithm
#' @param extending_ratio Minimum Ratio between windows with one different realisation to multiple in block-extending-algorithm
#' @param block_min_count_steps Number of steps till full filtering with block_min_count is done
#' @param snp_extending If TRUE use the SNP-extending-algorithm (step V) (default: TRUE)
#' @param max_extending_diff_snp Maximum number of SNPs with variants in SNP-extending-algorithm (step V; default: 0)
#' @param extending_ratio_snp Minimum ratio of SNPs with only one allele to those with variants (default: Inf)#' @param off_knot_addition If TRUE use off-variant-identification (default: FALSE)
#' @param off_knot_minimum_blocklength Minimum length of newly identified blocks (default: 10)
#' @param off_knot_minimum_blocksize Minimum number of individuals in newly identified blocks (default: 5)
#' @param raster Raster-width in the identification step (default: 5; recommended to be lower than off_knot_minimum_blocklength)
#' @param major_snp_calculation If TRUE calculate for major allele for each SNP in each block (default:TRUE)
#' @param off_knot_addition If TRUE identify additional blocks in regions not covered by the window cluster (default: FALSE)
#' @param prefilter If TRUE filter the dataset before detecting for block-structure (default: FALSE)
#' @param maf Minimum minor allel frequency in prefilter (default: 0.05)
#' @param equal_remove If TRUE filter out SNPs in perfect correlation to the next SNP (default: FALSE)
#' @param big_output If TRUE return not only blocklist but also blockinfo, dataset ect.
#' @param blockinfo_mode Structure of the groups in step I (default: 0 - Common haplos as major variants, 1- minimum number of groups)
#' @param c_dhm_mode If TRUE use high speed calculation with C in block_merging (default: TRUE)
#' @param intersect_func Selection of intersect function used in the computation (default: 2 - efficient C, 1- efficient R, 0- base::intersect)
#' @param fast_compiler If TRUE use compiler-package to enableJIT(3) (default: TRUE)
#' @param max_groups Maximum number of groups per window (adaptive window-size, default: 0 - use fixed window size)
#' @param recoding If TRUE change allele coding (Major allele "A", Minor allele "C")
#' @param recoding_notneeded Set to TRUE if dataset only contains "A", "C" to save comp. time
#' @param consider_multi If TRUE use multi-level edges to identify blocks (default: FALSE)
#' @param multi_min minimum number of haplotypes per multi transition to use in consider_multi (default: 5)
#' @param knoten_calc Set to 1 for an old version for calcuation of the data-list (NOT recommended)
#' @param blockinfo_mode2 Set to 4 for activation of NA-modelling
#' @param actual_snp_weight Set weight for difference between two alleles in a SNP (more than 1 possible base pair)
#' @param na_snp_weight Set weight for difference between NA and allele in a SNP (more than 1 possible base pair)
#' @param na_seq_weight Set weight for difference between NA and allele a loci with 1 possible base pair
#' @param weighting_length Weighting factor for length to determine major block (default: 1)
#' @param weighting_size Weighting factor for number of haplotypes in block to determine major block (default: 1)
#' @param target_coverage Target Coverage in the blocklist
#' @param max_iteration Maximum number of iterations to archive target coverage
#' @param min_step_size Minimum increase/decrase in min_block_size in target coverage fitting
#' @param target_stop Stop fitting target coverage fitting early
#' @param multi_window_mode Set TRUE so active the computation of multi window cluster with separate window_size/merging_error/mindestanteil
#' @param recalculate_biggest Set to FALSE to only calculate the number of major positions for those blocks that could be removed in each iteration (number of major blocks is only increasing when removing other blocks)
#' @export

block_calculation <- function(dhm, window_sequence=NULL, window_size=20, merging_error=1, minimum.blocksize=5, gap=10,
      mindestanteil=0.975, off_lines=5, mindestaehnlichkeit=0.99, merge_closeblock=FALSE,
      max_diff_l=1, max_diff_i=1, block_min_count=5000, bp_map=NULL, window_anchor_gens=NULL,
      consider_knoten=TRUE, consider_trans=TRUE, trans_min=5, subgroups=NULL, min_per_subgroup=0,
      consider_all=TRUE, save_allblock=TRUE, block_extending=TRUE,
      max_extending_diff=1, extending_ratio=20, block_min_count_steps=4,
      snp_extending=TRUE, max_extending_diff_snp=0, extending_ratio_snp=Inf,  major_snp_calculation=TRUE,
      off_knot_addition=FALSE, off_knot_minimum_blocklength=10, off_knot_minimum_blocksize=5,
      raster=5, at_least_one=TRUE,
      prefilter=FALSE, maf=0.00, equal_remove=FALSE,
      big_output=FALSE, blockinfo_mode=0, c_dhm_mode=TRUE,
      intersect_func=2, fast_compiler=TRUE,
      max_groups=0, recoding=FALSE, recoding_notneeded=FALSE,
      knoten_calc=2, consider_multi=FALSE, multi_min=5, blockinfo_mode2=3,
      actual_snp_weight = 5, na_snp_weight=2, na_seq_weight=0,
      weighting_length=1, weighting_size=1,
      recalculate_biggest=TRUE,
      target_coverage=NULL,
      max_iteration=5,
      min_step_size=25,
      target_stop=0.005,
      multi_window_mode=FALSE){

  if(multi_window_mode==TRUE){
    ncluster <- max(length(window_size), length(mindestanteil), merging_error)
    window_size <- rep(window_size,length.out=ncluster)
    mindestanteil <- rep(mindestanteil, length.out=ncluster)
    merging_error <- rep(merging_error, length.out=ncluster)
  } else{
    ncluster <- 1
  }

  if(merge_closeblock==TRUE && length(unique(window_size))){
    merge_closeblock <- FALSE
    print("Closeblock-Merging only for single window size")
  }

  window_sequence_list <- list()

  if(length(window_sequence)!=ncluster){
    for(index in 1:ncluster){
      window_sequence_list[[index]] <- window_sequence
    }
  } else{
    window_sequence_list <- window_sequence
  }
  window_sequence_list[[ncluster+1]] <- "placeholder"



  if(fast_compiler){
    requireNamespace("compiler")
    enableJIT(3)

  }

  if(intersect_func==0){
    intersect_func <- base::intersect
  } else if(intersect_func==1){
    intersect_func <- intersect_own
  } else if(intersect_func==2){
    intersect_func <- HaploBlocker::intersect
  } else if(intersect_func==3){
    intersect_func <- intersect_own3
  }

  if(is.data.frame(dhm)){
    dhm <- as.matrix(dhm)
  }

  if(prefilter==TRUE){
    dhm <- dataset_filter(dhm, maf, equal_remove)
  }

  if(recoding==TRUE){
    if(recoding_notneeded==TRUE){
      recoding <- TRUE
    } else{
      for(index in 1:nrow(dhm)){
        check1 <- dhm[index,]==dhm[index,1]
        dhm[index, check1] <- "A"
        dhm[index, -(check1)*1:indi] <- "C"
      }
    }
  }

  # Martins Erweiterung
  unique.dhm <- unique(as.vector(dhm))
  fixcoding(unique.dhm)
  c_dhm <- codeSNPs(dhm)

  blockinfo <- list()
  data <- list()
  dataset <- list()
  partial_blocklist <- list()
  nwindow <- rep(0, ncluster)

  for(index in 1:ncluster){
    if(blockinfo_mode2==4){
      blockinfo_out <- blockinfo_calculation4(dhm, window_sequence= window_sequence_list[[index]], window_anchor_gens = window_anchor_gens, blockinfo_mode=blockinfo_mode,
                                              window_size = window_size[index], merging_error = merging_error[index], bp_map = bp_map, at_least_one=at_least_one,
                                              actual_snp_weight=actual_snp_weight, na_snp_weight=na_snp_weight, na_seq_weight= na_seq_weight)

    } else{
      blockinfo_out <- blockinfo_calculation3(dhm, window_sequence= window_sequence_list[[index]], window_anchor_gens = window_anchor_gens, blockinfo_mode=blockinfo_mode,
                                              window_size = window_size[index], merging_error = merging_error[index], bp_map = bp_map, at_least_one=at_least_one,
                                              c_dhm=c_dhm, max_groups=max_groups)

    }

    blockinfo[[index]] <- blockinfo_out[[1]]
    window_sequence_list[[index]] <- blockinfo_out[[2]]


    if(max_groups>0){
      print(paste("Generated:", nrow(window_sequence_list[[index]]), "windows"))
      print(paste("With size: Max", max(window_sequence_list[[index]][,3]), "Min", min(window_sequence_list[[index]][,3]), "Avg", round(mean(window_sequence_list[[index]][,3])*100)/100))
    }


    if(knoten_calc==1){
      data[[index]] <- knoten_calculation(blockinfo[[index]], window_sequence_list[[index]])
    } else{
      data[[index]] <- knoten_calculation2(blockinfo[[index]], window_sequence_list[[index]])

    }

    data[[index]] <- simple_merge(data[[index]], intersect_func=intersect_func)
    data[[index]] <- calculate_transition(data[[index]], intersect_func=intersect_func)


      nwindow[index] <- nrow(window_sequence_list[[index]])
      indi <- sum(blockinfo[[index]][[1]][[1]])
      knoten <- length(data[[index]])
      iteration <- 1
      print("Start_CrossMerging_full")
      a <- start_end_block(data[[index]])
      a_old <- NULL
      while(length(a_old)==0 || (nrow(a_old)!=nrow(a)) || prod(a_old==a)==0){
        a_old <- a
        print(paste("Iteration", iteration, ":", knoten, "Knoten"))
        data[[index]] <- crossmergev2(data[[index]], indi, nwindow[index], a, intersect_func=intersect_func)
        data[[index]] <- simple_merge_probv2(data[[index]], indi, nwindow[index])
        a <- start_end_block(data[[index]])
        knoten <- length(data[[index]])
        iteration <- iteration + 1
      }


      knoten <- length(data[[index]])
      iteration <- 1
      print("Start_IgnoreSmall")
      while(iteration==1 || (nrow(a_old)!=nrow(a)) || prod(a_old==a)==0){
        a_old <- a
        print(paste("Iteration", iteration, ":", knoten, "Knoten"))
        data[[index]] <- ignore_small_knotsv2(data[[index]], indi, nwindow[index], minimum.blocksize, gap, intersect_func=intersect_func)
        data[[index]] <- simple_merge_probv2(data[[index]], indi, nwindow[index], intersect_func=intersect_func)
        data[[index]] <- crossmergev2(data[[index]], indi, nwindow[index], intersect_func=intersect_func)
        data[[index]] <- simple_merge_probv2(data[[index]], indi, nwindow[index], intersect_func=intersect_func)
        a<- start_end_block(data[[index]])
        knoten <- length(data[[index]])
        iteration <- iteration + 1
      }

    dataset[[index]] <- block_dataset_construction(blockinfo[[index]], indi, nwindow[index])

    partial_blocklist[[index]] <- identify_blocks3(data[[index]], indi, nwindow[index], mindestanteil[index], trans_min=trans_min, subgroups=subgroups,
                                                 consider_knoten=consider_knoten, consider_trans=consider_trans, min_per_subgroup=min_per_subgroup,
                                                 intersect_func=intersect_func, consider_multi=consider_multi, multi_min=multi_min)
  }

  blocklist <- list()
  for(index in 1:ncluster){
    for(index2 in 1:length(partial_blocklist[[index]])){
      blocklist[[length(blocklist)+1]] <- partial_blocklist[[index]][[index2]]
      blocklist[[length(blocklist)]][[12]] <- index
    }
  }

  # STRINGENZ ZU DATASET (ncol/nrow)


  current_iteration <- max_iteration
  if(length(target_coverage)>0){
    current_iteration <- 1
    blocklist_start <- blocklist
    block_min_counts <- numeric(max_iteration)
    coverage_results <- numeric(max_iteration)
    block_min_counts[1] <- block_min_count
    blocklists <- list()

  }
  stop_iteration <- FALSE
  while(current_iteration <= max_iteration && !stop_iteration){

    if(length(target_coverage)>0){
      block_min_count <- block_min_counts[current_iteration]
      blocklist <- blocklist_start
    }

    knoten <- length(blocklist)
    iteration <- 1
    print("Start_Blockmerging")
    helper <- blocklist_startend(blocklist, type="snp")
    helper_old <- NULL
    while(length(helper_old)==0 || (nrow(helper_old)!=nrow(helper)) || prod(helper_old==helper)==0){
      print(paste("Iteration", iteration, ":", knoten, "Bloecke"))
      helper_old <- helper
      blocklist <- block_merging(blocklist, blockinfo, dataset, dhm, indi, nwindow, window_sequence_list, off_lines, mindestaehnlichkeit=mindestaehnlichkeit,
                                 consider_all=consider_all, minimum.blocksize=minimum.blocksize, save_allblock=save_allblock, helper=helper,
                                 c_dhm=c_dhm, c_dhm_mode=c_dhm_mode, intersect_func=intersect_func,
                                 min_per_subgroup=min_per_subgroup)

      if(merge_closeblock==TRUE){
        blocklist <- block_closeblock_merging(blocklist, blockinfo, indi, nwindow, max_diff_l, max_diff_i, intersect_func=intersect_func,
                                              dataset=dataset)
      }
      if(block_min_count>(0)){
        if(block_min_count_steps>1){
          blocklist <- blockinfo_biggest(blocklist, nwindow, indi, block_min_count=min(iteration-1, block_min_count_steps-1)/(block_min_count_steps-1)*block_min_count,
                                         weighting_length=weighting_length, weighting_size=weighting_size,
                                         recalculate_biggest=recalculate_biggest, window_size=window_size)
        } else{
          blocklist <- blockinfo_biggest(blocklist, nwindow, indi, block_min_count=block_min_count, weighting_length=weighting_length, weighting_size=weighting_size,
                                         recalculate_biggest=recalculate_biggest, window_size=window_size)
        }
      }
      helper <- blocklist_startend(blocklist, type="snp")
      knoten <- length(blocklist)
      iteration <- iteration + 1
    }

    if(block_extending==TRUE){
      knoten <- length(blocklist)
      iteration <- 1
      extentions_done <- 0
      print("Start_Blockextending")
      while(iteration==1 || (nrow(helper_old)!=nrow(helper)) || prod(helper_old==helper)==0 || extentions_done > 0){
        helper_old <- helper
        print(paste("Iteration", iteration, ":", knoten, "Bloecke; ", extentions_done, "Blockerweiterungen"))
        blocklist_out <- extend_block(blocklist, indi, nwindow, max_extending_diff=max_extending_diff,
                                      extending_ratio=extending_ratio, dataset=dataset, window_sequence_list=window_sequence_list)
        blocklist <- blocklist_out[[1]]
        extentions_done <- blocklist_out[[2]]

        blocklist <- block_merging(blocklist, blockinfo, dataset, dhm, indi, nwindow, window_sequence_list, off_lines, mindestaehnlichkeit=mindestaehnlichkeit,
                                   consider_all=consider_all, minimum.blocksize=minimum.blocksize, save_allblock=save_allblock, c_dhm=c_dhm,
                                   c_dhm_mode=c_dhm_mode , intersect_func=intersect_func,
                                   min_per_subgroup=min_per_subgroup)

        if(merge_closeblock==TRUE){
          blocklist <- block_closeblock_merging(blocklist, blockinfo, indi, nwindow, max_diff_l, max_diff_i, intersect_func=intersect_func,
                                                dataset=dataset)
        }
        if(block_min_count_steps>1){
          blocklist <- blockinfo_biggest(blocklist, nwindow, indi, block_min_count=min(iteration-1, block_min_count_steps-1)/(block_min_count_steps-1)*block_min_count,
                                         weighting_length=weighting_length, weighting_size=weighting_size,
                                         recalculate_biggest=recalculate_biggest, window_size=window_size)
        } else{
          blocklist <- blockinfo_biggest(blocklist, nwindow, indi, block_min_count=block_min_count,
                                         weighting_length=weighting_length, weighting_size=weighting_size,
                                         recalculate_biggest=recalculate_biggest, window_size=window_size)
        }

        helper <- blocklist_startend(blocklist, type="snp")
        knoten <- length(blocklist)
        iteration <- iteration + 1
      }
    }

    if(snp_extending==TRUE){
      # Fuege anfang und ende gleiche SNPs im Randbereich hinzu/ab
      blocklist <- extend_snp(blocklist, indi, nwindow, dhm, window_sequence_list, bp_map = bp_map,
                              max_extending_diff_snp=max_extending_diff_snp, extending_ratio_snp=extending_ratio_snp)

    }



    if(off_knot_addition==TRUE){
      print("OFF_Knot_addition Only for single nwindow")
      blocklist <- add_offknot(blocklist=blocklist, dataset[[1]], indi, nwindow[[1]], window_sequence=window_sequence_list[[1]], bp_map=bp_map,
                               off_knot_minimum_blocklength=off_knot_minimum_blocklength,
                               off_knot_minimum_blocksize=off_knot_minimum_blocksize, raster=raster)
      knoten <- length(blocklist)
      iteration <- 1
      print("Start_Blockmerging")
      helper <- blocklist_startend(blocklist, type="snp")
      helper_old <- NULL
      while(length(helper_old)==0 || (nrow(helper_old)!=nrow(helper)) || prod(helper_old==helper)==0){
        print(paste("Iteration", iteration, ":", knoten, "Bloecke"))
        helper_old <- helper
        blocklist <- block_merging(blocklist, blockinfo, dataset, dhm, indi, nwindow, window_sequence_list, off_lines, mindestaehnlichkeit=mindestaehnlichkeit,
                                   consider_all=consider_all, minimum.blocksize=minimum.blocksize, save_allblock=save_allblock, helper=helper,
                                   c_dhm=c_dhm, c_dhm_mode=c_dhm_mode, intersect_func=intersect_func,
                                   min_per_subgroup=min_per_subgroup)

        if(merge_closeblock==TRUE){
          blocklist <- block_closeblock_merging(blocklist, blockinfo, indi, nwindow, max_diff_l, max_diff_i, intersect_func=intersect_func,
                                                dataset=dataset)
        }
        if(block_min_count>(0)){
          if(block_min_count_steps>1){
            blocklist <- blockinfo_biggest(blocklist, nwindow, indi, block_min_count=min(iteration-1, block_min_count_steps-1)/(block_min_count_steps-1)*block_min_count,
                                           weighting_length=weighting_length, weighting_size=weighting_size,
                                           recalculate_biggest=recalculate_biggest, window_size=window_size)
          } else{
            blocklist <- blockinfo_biggest(blocklist, nwindow, indi, block_min_count=block_min_count, weighting_length=weighting_length, weighting_size=weighting_size,
                                           recalculate_biggest=recalculate_biggest, window_size=window_size)
          }
        }
        helper <- blocklist_startend(blocklist, type="snp")
        knoten <- length(blocklist)
        iteration <- iteration + 1
      }

      if(block_extending==TRUE){

        knoten <- length(blocklist)
        iteration <- 1
        extentions_done <- 0
        print("Start_Blockextending")
        while(iteration==1 || (nrow(helper_old)!=nrow(helper)) || prod(helper_old==helper)==0 || extentions_done > 0){
          helper_old <- helper
          print(paste("Iteration", iteration, ":", knoten, "Bloecke; ", extentions_done, "Blockerweiterungen"))
          blocklist_out <- extend_block(blocklist, indi, nwindow, max_extending_diff=max_extending_diff,
                                        extending_ratio=extending_ratio, dataset=dataset, window_sequence_list=window_sequence_list)
          blocklist <- blocklist_out[[1]]
          extentions_done <- blocklist_out[[2]]

          blocklist <- block_merging(blocklist, blockinfo, dataset, dhm, indi, nwindow, window_sequence_list, off_lines, mindestaehnlichkeit=mindestaehnlichkeit,
                                     consider_all=consider_all, minimum.blocksize=minimum.blocksize, save_allblock=save_allblock, c_dhm=c_dhm,
                                     c_dhm_mode=c_dhm_mode , intersect_func=intersect_func,
                                     min_per_subgroup=min_per_subgroup)

          if(merge_closeblock==TRUE){
            blocklist <- block_closeblock_merging(blocklist, blockinfo, indi, nwindow, max_diff_l, max_diff_i, intersect_func=intersect_func,
                                                  dataset=dataset)
          }
          if(block_min_count_steps>1){
            blocklist <- blockinfo_biggest(blocklist, nwindow, indi, block_min_count=min(iteration-1, block_min_count_steps-1)/(block_min_count_steps-1)*block_min_count,
                                           weighting_length=weighting_length, weighting_size=weighting_size,
                                           recalculate_biggest=recalculate_biggest, window_size=window_size)
          } else{
            blocklist <- blockinfo_biggest(blocklist, nwindow, indi, block_min_count=block_min_count,
                                           weighting_length=weighting_length, weighting_size=weighting_size,
                                           recalculate_biggest=recalculate_biggest, window_size=window_size)
          }

          helper <- blocklist_startend(blocklist, type="snp")
          knoten <- length(blocklist)
          iteration <- iteration + 1
        }
      }

      if(snp_extending==TRUE){
        # Fuege anfang und ende gleiche SNPs im Randbereich hinzu/ab
        blocklist <- extend_snp(blocklist, indi, nwindow, dhm, window_sequence_list, bp_map = bp_map,
                                max_extending_diff_snp=max_extending_diff_snp, extending_ratio_snp=extending_ratio_snp)

      }


    }

    if(length(target_coverage)>0){

      blocklists[[current_iteration]] <- blocklist
      t <- coverage_test(blocklist, indi, type="snp")
      coverage_results[current_iteration] <- mean(t)
      prev_cov <- coverage_results[coverage_results>0]
      prev_block <- block_min_counts[block_min_counts>0]
      if(min(prev_cov) > target_coverage){
        block_min_counts[current_iteration + 1] <- max(max(prev_block)*2, max(prev_block)+min_step_size)
      } else if(max(prev_cov) < target_coverage){
        block_min_counts[current_iteration + 1] <- max(min(ceiling(min(prev_block)/2), min(prev_block) - min_step_size),1)
      } else{
        if(current_iteration==1){
          stop_iteration <- TRUE
        } else{
          pos <- ((prev_cov-target_coverage)>0) * (1:length(prev_cov))
          neg <- ((prev_cov-target_coverage)<0) * (1:length(prev_cov))
          ordering <- sort(abs(prev_cov-target_coverage), index.return=TRUE)$ix
          neg <- base::intersect(ordering, neg)
          pos <- base::intersect(ordering, pos)
          block_min_counts[current_iteration + 1] <- ceiling(mean(prev_block[c(pos[1], neg[1])]))
        }

      }
      print(paste("Finish Target_Coverage Iteration", current_iteration))
      print(paste("Used block_min_count:", block_min_counts[current_iteration]))
      print(paste("Achieved Coverage:", coverage_results[current_iteration]))
      if(current_iteration<max_iteration){
        print(paste("Start next iteration using block_min_count:", block_min_counts[current_iteration+1]))
      }

      if(abs(coverage_results[current_iteration] - target_coverage)<target_stop){
        stop_iteration <- TRUE
      }




    }
    current_iteration <- current_iteration + 1

  }

  if(length(target_coverage)>0){
    print(paste("Final Iteration using ",  block_min_counts[current_iteration-1]))
    print(paste("Achieved Coverage:", coverage_results[current_iteration-1]))
    take <- which.max(-abs(coverage_results-target_coverage))[1]
    blocklist <- blocklists[[take]]
  }


  if(major_snp_calculation==TRUE){
    blocklist <- major_snp_calculation(blocklist, dhm, recoding=recoding)
  }



  # Leere Blöcke - wesentlich kleinere Blöcke nicht direkt entfernen sondern verkleinerte versionen behalten.
  if(big_output==FALSE){
    return(blocklist)
  } else{
    return(list(blocklist, dataset, data, blockinfo, indi, nwindow))
  }

}
