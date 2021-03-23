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
#' @param node_min minimum number of haplotypes per block (default: 5)
#' @param gap remove haplotypes in nodes adjacent to nodes with less than min_node haplotypes in it (default: 10 windows)
#' @param min_share minimum percentage of transition to the same block for extension (default: 0.95, step III)
#' @param off_lines minimum number of haplotypes to looose in the creation of a bigger block (default: 5)
#' @param min_similarity minimum rate of the same SNPs to be added to the block (default: 0.99)
#' @param merge_closeblock If TRUE merge adjecting block with similar haplotypes (default: FALSE)
#' @param max_diff_l maximum number of windows with different haplotypes inbetween (default: 1)
#' @param max_diff_i maximum number of individuals in only one of the two blocks (default: 1)
#' @param min_majorblock minimum of cells in the dataset a block is the biggest covering (default: 5000)
#' @param consider_nodes Use nodes to identify blocks (default: TRUE)
#' @param consider_edge Use edges between nodes to identify blocks (default: TRUE)
#' @param edge_min minimum number of haplotypes per transition to use in consider_edge (default: 5)
#' @param subgroups possible subgroups to consider in the block identification process (default: NULL - list(1:indi))
#' @param min_per_subgroup minimum number of haplotypes per block per subgroup (default: 0)
#' @param subgroup_exception allow for not all subgroups to include a block
#' @param consider_all If TRUE always haplotypes which are not in the node to be in a generated block
#' @param save_allblock If TRUE keep all haplotypes with all windows according to a block (even under min_similarity)
#' @param block_extending If TRUE use the window-extending algorithm (step V)
#' @param max_extending_diff Maximum number of windows with different realisation in the block-extending-algorithm
#' @param extending_ratio Minimum Ratio between windows with one different realisation to multiple in block-extending-algorithm
#' @param min_majorblock_steps Number of steps till full filtering with min_majorblock is done
#' @param snp_extending If TRUE use the SNP-extending-algorithm (step V) (default: TRUE)
#' @param max_extending_diff_snp Maximum number of SNPs with variants in SNP-extending-algorithm (step V; default: 0)
#' @param extending_ratio_snp Minimum ratio of SNPs with only one allele to those with variants (default: Inf)#'
#' @param off_node_addition If TRUE identify additional blocks in regions not covered by the window cluster (default: FALSE)
#' @param off_node_minimum_blocklength Minimum length of newly identified blocks (default: 10)
#' @param off_node_minimum_blocksize Minimum number of individuals in newly identified blocks (default: 5)
#' @param raster Raster-width in the identification step (default: 5; recommended to be lower than off_node_minimum_blocklength)
#' @param major_snp_calculation If TRUE calculate for major allele for each SNP in each block (default:TRUE)
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
#' @param blockinfo_mode_na Set TRUE for activation of NA-modelling
#' @param actual_snp_weight Set weight for difference between two alleles in a SNP (more than 1 possible base pair)
#' @param na_snp_weight Set weight for difference between NA and allele in a SNP (more than 1 possible base pair)
#' @param na_seq_weight Set weight for difference between NA and allele a loci with 1 possible base pair
#' @param weighting_length Weighting factor for length to determine major block (default: 1)
#' @param weighting_size Weighting factor for number of haplotypes in block to determine major block (default: 1)
#' @param target_coverage Target Coverage in the blocklist
#' @param max_iteration Maximum number of iterations to archive target coverage
#' @param min_step_size Minimum increase/decrase in min_block_size in target coverage fitting
#' @param target_stop Stop fitting target coverage fitting early
#' @param multi_window_mode Set TRUE so active the computation of multi window cluster with separate window_size/merging_error/min_share
#' @param recalculate_biggest Set to FALSE to only calculate the number of major positions for those blocks that could be removed in each iteration (number of major blocks is only increasing when removing other blocks)
#' @param adaptive_mode Set to TRUE to use the predefined adoptive mode ( window_size = c(5,10,20,50), target_coverage=0.9)
#' @param developer_mode Set to TRUE to not delete internal computational stuff in the blocklist
#' @param parallel_window Maximum number of markers to process jointly (default: Inf)
#' @param window_overlap Overlap between windows in parallel mode (default: 0)
#' @param window_cores Number of cores to use in parallel mode (default: 1)
#' @param double_share Extended-block-identification share t
#' @param min_reduction_cross Minimum number of nodes reduction to continue SM, SG cycle (default: -Inf, till fixation)
#' @param min_reduction_neglet Minimum number of nodes reduction to continue NN, SM, SG, SM cycle (default: -Inf, till fixation)
#' @param early_remove Remove Nodes with little number of haplotypes before SM, SG cycle (default: FALSE)
#' @param node_min_early Minimum number of haplotypes per node before SM, SG cycle (default: NULL)
#' @param overlap_remove If set to TRUE the obtained Haplotype Library will have no overlapping blocks.
#' @param deletion_count If TRUE 0s are handles as deletions. Not increasing rating // counted as major positions
#' @param na_value Number/Character variable that is coding NA (default: -9)
#' @param verbose Set to FALSE to not display any prints
#' @param inbred Set to TRUE when working with inbred material and importing genomic data from file
#' @examples
#' data(ex_maze)
#' blocklist <- block_calculation(ex_maze)
#' @export
#' @return haplotype library
#'


block_calculation <- function(dhm, window_sequence=NULL, window_size=20, merging_error=1, node_min=5, gap=10,
                              min_share=0.975, off_lines=2, min_similarity=0.99, merge_closeblock=FALSE,
                              max_diff_l=1, max_diff_i=1, min_majorblock=5000, bp_map=NULL, window_anchor_gens=NULL,
                              consider_nodes=TRUE, consider_edge=TRUE, edge_min=5, subgroups=NULL, min_per_subgroup=0, subgroup_exception=0,
                              consider_all=TRUE, save_allblock=TRUE, block_extending=TRUE,
                              max_extending_diff=1, extending_ratio=20, min_majorblock_steps=4,
                              snp_extending=TRUE, max_extending_diff_snp=0, extending_ratio_snp=Inf,  major_snp_calculation=TRUE,
                              off_node_addition=FALSE, off_node_minimum_blocklength=10, off_node_minimum_blocksize=5,
                              raster=5, at_least_one=TRUE,
                              prefilter=FALSE, maf=0.00, equal_remove=FALSE,
                              big_output=FALSE, blockinfo_mode=0, c_dhm_mode=TRUE,
                              intersect_func=TRUE, fast_compiler=TRUE,
                              max_groups=0, recoding=FALSE, recoding_notneeded=FALSE,
                              consider_multi=FALSE, multi_min=5, blockinfo_mode_na=FALSE,
                              actual_snp_weight = 5, na_snp_weight=2, na_seq_weight=0,
                              weighting_length=1, weighting_size=1,
                              recalculate_biggest=TRUE,
                              target_coverage=NULL,
                              max_iteration=10,
                              min_step_size=25,
                              target_stop=0.001,
                              multi_window_mode=FALSE,
                              adaptive_mode=FALSE,
                              developer_mode=FALSE,
                              double_share=1,
                              early_remove=FALSE,
                              node_min_early=NULL,
                              min_reduction_cross=-Inf,
                              min_reduction_neglet=-Inf,
                              parallel_window=Inf,
                              window_overlap=0,
                              window_cores=1, overlap_remove=FALSE,
                              na_value=(-9),
                              deletion_count=FALSE,
                              verbose=TRUE,
                              inbred = FALSE){

  if(adaptive_mode==TRUE){
    multi_window_mode <- TRUE
    if((length(window_size)==1 && window_size==20)){
      window_size <- c(5,10,20,50)
      if(overlap_remove){
        window_size <- 20
        warning("Adaptive mode only use a single window size (20) when using overlap removal.")
      }
    }
    if(length(target_coverage)==0){
      target_coverage <- 0.9
    }
  }

  if(length(window_size)>1 && overlap_remove){
    if(adaptive_mode){
      warning("Adaptive mode will use multiple different window sizes")
    }
    stop("The use of multiple window sizes and overlap removal at the same time is not supported!")
  }

  if(length(subgroups)>0){
    for(index in 1:length(subgroups)){
      subgroups[[index]] <- sort(subgroups[[index]])
    }
  }


  {
    # foreach variables
    indexb <- NULL
    if (requireNamespace("foreach", quietly = TRUE)) {
      `%dopar%` <- foreach::`%dopar%`
    }

  }
  if(multi_window_mode==FALSE && (length(window_size)>1 || length(merging_error)>1 || length(min_share)>1)){
    if(verbose) cat("Active multi_window_mode\n")
    multi_window_mode <- TRUE
  }


  if(multi_window_mode){
    ncluster <- max(length(window_size), length(min_share), merging_error)
    window_size <- rep(window_size,length.out=ncluster)
    min_share <- rep(min_share, length.out=ncluster)
    merging_error <- rep(merging_error, length.out=ncluster)
  } else{
    ncluster <- 1
  }


  if(length(dhm)==1){

    data_file <- data_import(dhm, inbred=inbred, verbose=verbose)
    dhm <- data_file[[1]]
    if(length(bp_map)==0){
      bp_map <- data_file[[2]]
    }
    rm(data_file)
    if(verbose) cat("Data import successful. \n")
  }

  if(length(bp_map)>0){
    bp_map <- as.numeric(bp_map)
  } else{
    bp_map <- rep(0, nrow(dhm))
  }

  if(length(window_size)==1 && nrow(dhm) && window_size >= nrow(dhm)){
    warning("Your dataset is extremely small - resulting in just one window for the window cluster!")
    warning("window_size has been set to 1")
    window_size <- 1

  }

  if(nrow(dhm)>parallel_window){
    dhm_list <- list()
    snp_start <- NULL
    bp_map_list <- list()
    n_windows <- max(1,ceiling((nrow(dhm)-window_overlap)/parallel_window))
    for(index in 1:n_windows){
      snp_start <- c(snp_start, (index-1)*parallel_window)
      if(index!=n_windows){
        dhm_list[[index]] <- dhm[1:(parallel_window+window_overlap) + (index-1)*parallel_window,]
        if(length(bp_map>0)){
          bp_map_list[[index]] <- bp_map[1:(parallel_window+window_overlap) + (index-1)*parallel_window]
        } else{
          bp_map_list[[index]] <- NULL
        }
      } else{
        dhm_list[[index]] <- dhm[(1+(index-1)*parallel_window):nrow(dhm),]
        if(length(bp_map>0)){
          bp_map_list[[index]] <- bp_map[(1+(index-1)*parallel_window):nrow(dhm)]
        } else{
          bp_map_list[[index]] <- NULL
        }
      }

    }
    bp_map_list[[index+1]] <- "placeholder"
    blocklist_list <- NULL

    block_calculation_par <- function(dhm){
      blockl <- block_calculation(dhm,
                        window_sequence=window_sequence, window_size=window_size,
                        merging_error=merging_error, node_min=node_min, gap=gap,
                        min_share=min_share, off_lines=off_lines, min_similarity=min_similarity,
                        merge_closeblock=merge_closeblock,
                        max_diff_l=max_diff_l, max_diff_i=max_diff_i, min_majorblock=min_majorblock,
                        bp_map=NULL, window_anchor_gens=window_anchor_gens,
                        consider_nodes=consider_nodes, consider_edge=consider_edge,
                        edge_min=edge_min, subgroups=subgroups, min_per_subgroup=min_per_subgroup,
                        subgroup_exception=subgroup_exception,
                        consider_all=consider_all, save_allblock=save_allblock,
                        block_extending=block_extending,
                        max_extending_diff=max_extending_diff, extending_ratio=extending_ratio,
                        min_majorblock_steps=min_majorblock_steps,
                        snp_extending=snp_extending, max_extending_diff_snp=max_extending_diff_snp,
                        extending_ratio_snp=extending_ratio_snp,  major_snp_calculation=major_snp_calculation,
                        off_node_addition=off_node_addition,
                        off_node_minimum_blocklength=off_node_minimum_blocklength,
                        off_node_minimum_blocksize=off_node_minimum_blocksize,
                        raster=raster, at_least_one=at_least_one,
                        prefilter=prefilter, maf=maf, equal_remove=equal_remove,
                        big_output=FALSE, blockinfo_mode=blockinfo_mode, c_dhm_mode=c_dhm_mode,
                        intersect_func=intersect_func, fast_compiler=fast_compiler,
                        max_groups=max_groups, recoding=recoding, recoding_notneeded=recoding_notneeded,
                        consider_multi=consider_multi, multi_min=multi_min,
                        blockinfo_mode_na=blockinfo_mode_na,
                        actual_snp_weight = actual_snp_weight, na_snp_weight=na_snp_weight,
                        na_seq_weight=na_seq_weight,
                        weighting_length=weighting_length, weighting_size=weighting_size,
                        recalculate_biggest=recalculate_biggest,
                        target_coverage=target_coverage,
                        max_iteration=max_iteration,
                        min_step_size=min_step_size,
                        target_stop=target_stop,
                        multi_window_mode=multi_window_mode,
                        adaptive_mode=adaptive_mode,
                        developer_mode=TRUE,
                        double_share=double_share,
                        early_remove=early_remove,
                        node_min_early=node_min_early,
                        min_reduction_cross=min_reduction_cross,
                        min_reduction_neglet=min_reduction_neglet,
                        deletion_count = deletion_count,
                        na_value = na_value,
                        overlap_remove = overlap_remove,
                        verbose = verbose,
                        inbred = inbred

                        )
      return(blockl)
    }

    if(Sys.info()[['sysname']]=="Windows"){
      doParallel::registerDoParallel(cores=window_cores)
      blocklist_list <- foreach::foreach(indexb=1:length(dhm_list), .packages = "HaploBlocker") %dopar% {
        element <- block_calculation_par(dhm_list[[indexb]])
        element
      }
      doParallel::stopImplicitCluster()
    } else{
      blocklist_list <- parallel::mclapply(dhm_list, block_calculation_par, mc.cores=window_cores)
    }

   rm(dhm_list)


    # modify Start/End
    for(index in 2:length(blocklist_list)){
      for(index2 in 1:length(blocklist_list[[index]])){
        if(length(bp_map_list[[index]])>0){
          blocklist_list[[index]][[index2]][[2]]$bp <- bp_map_list[[index]][blocklist_list[[index]][[index2]][[2]]$snp]
          blocklist_list[[index]][[index2]][[3]]$bp <- bp_map_list[[index]][blocklist_list[[index]][[index2]][[3]]$snp]
        }
        blocklist_list[[index]][[index2]][[2]]$snp <- blocklist_list[[index]][[index2]][[2]]$snp + snp_start[index]
        blocklist_list[[index]][[index2]][[3]]$snp <- blocklist_list[[index]][[index2]][[3]]$snp + snp_start[index]
        blocklist_list[[index]][[index2]][[2]]$window <- blocklist_list[[index]][[index2]][[2]]$window + snp_start[index]/window_size[blocklist_list[[index]][[index2]][[12]]]
        blocklist_list[[index]][[index2]][[3]]$window <- blocklist_list[[index]][[index2]][[3]]$window + snp_start[index]/window_size[blocklist_list[[index]][[index2]][[12]]]
      }
    }
    blocklist <- NULL
    for(index in 1:length(blocklist_list)){
      blocklist <- c(blocklist, blocklist_list[[index]])
    }


    ### PRESENT DATA NEEDS TO BE RECALCULATED HERE!
    nwindow <- nrow(dhm) / window_size
    present_data <- matrix(0, nrow= ncol(dhm), ncol = nwindow)

    if(deletion_count){
      dhm_p <- dhm!=na_value
      while(nrow(dhm)< window_size * nwindow){
        dhm_p <- rbind((dhm!= na_value),0)
      }

      for(temp1 in 1:window_size){
        present_data <- present_data + t(dhm_p[1:nwindow*window_size-temp1+1,])
      }
      present_data <- present_data / window_size
    }


    if(window_overlap>0){
      blocklist <- blocklist_reorder(blocklist, node_min)
      blocklist <- blockinfo_biggest(blocklist, min_majorblock=min_majorblock, weighting_length=weighting_length, weighting_size=weighting_size,
                                     recalculate_biggest=recalculate_biggest, window_size=window_size,
                                     deletion_count=deletion_count, present_data=present_data)

      if(overlap_remove){
        warning("Overlapping segments between windows can contain overlapping blocks. If that is not wanted set window_overlap = 0 or not run in parallel!")
      }

    }

    if(big_output && developer_mode){
      return(list(blocklist, NULL, NULL, NULL, NULL, NULL))
    } else if(big_output){
      return(list(blocklist, NULL, NULL, NULL, NULL, NULL, NULL, present_data))
    } else{
      return(blocklist)
    }

  }

  if(merge_closeblock==TRUE && length(unique(window_size))){
    merge_closeblock <- FALSE
    if(verbose) cat("Closeblock-Merging only for single window size\n")
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

  if(is.logical(intersect_func) && intersect_func){
    intersect_func <- HaploBlocker::intersect_c
  } else if(is.logical(intersect_func) && !intersect_func){
    intersect_func <- base::intersect
  } else{
    intersect_func <- intersect_func
  }

  if(is.data.frame(dhm)){
    dhm <- as.matrix(dhm)
  }
  if(sum(is.na(dhm))>0){
    dhm[is.na(dhm)] <- na_value
  }

  if(prefilter==TRUE){
    dhm_1 <- dataset_filter(dhm, maf, equal_remove, bp_map=bp_map)
    dhm <- dhm_1[[1]]
    bp_map <- dhm_1[[2]]
  }

  if(recoding==TRUE){
    if(recoding_notneeded==TRUE){
      recoding <- TRUE
    } else{
      for(index in 1:nrow(dhm)){
        check1 <- dhm[index,]==dhm[index,1]
        dhm[index, check1] <- "A"
        dhm[index, -(check1)*1:ncol(dhm)] <- "C"
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
    if(blockinfo_mode_na){
      blockinfo_out <- blockinfo_calculation_na(dhm, window_sequence= window_sequence_list[[index]], window_anchor_gens = window_anchor_gens, blockinfo_mode=blockinfo_mode,
                                                window_size = window_size[index], merging_error = merging_error[index], bp_map = bp_map, at_least_one=at_least_one,
                                                actual_snp_weight=actual_snp_weight, na_snp_weight=na_snp_weight, na_seq_weight= na_seq_weight,
                                                verbose=verbose)

    } else{
      blockinfo_out <- blockinfo_calculation(dhm, window_sequence= window_sequence_list[[index]], window_anchor_gens = window_anchor_gens, blockinfo_mode=blockinfo_mode,
                                             window_size = window_size[index], merging_error = merging_error[index], bp_map = bp_map, at_least_one=at_least_one,
                                             c_dhm=c_dhm, max_groups=max_groups,
                                             verbose=verbose)

    }

    blockinfo[[index]] <- blockinfo_out[[1]]
    window_sequence_list[[index]] <- blockinfo_out[[2]]


    if(max_groups>0){
      if(verbose) cat(paste("Generated:", nrow(window_sequence_list[[index]]), "windows\n"))
      if(verbose) cat(paste("With size: Max", max(window_sequence_list[[index]][,3]), "Min", min(window_sequence_list[[index]][,3]), "Avg", round(mean(window_sequence_list[[index]][,3])*100)/100),"\n")
    }

    data[[index]] <- nodes_calculation(blockinfo[[index]], window_sequence_list[[index]], verbose = verbose)

    data[[index]] <- simple_merge(data[[index]], intersect_func=intersect_func, verbose = verbose)
    data[[index]] <- calculate_transition(data[[index]], intersect_func=intersect_func)



    nwindow[index] <- nrow(window_sequence_list[[index]])
    indi <- sum(blockinfo[[index]][[1]][[1]])

    nodes <- length(data[[index]])

    if(early_remove){
      if(verbose) cat("Start_Early_remove\n")
      if(verbose) cat(paste("Starting-Nodes:", nodes,"\n"))
      if(length(node_min_early)==0){
        node_min_early <- node_min
      }
      data[[index]] <- ignore_small_nodes(data[[index]], indi, nwindow[index], node_min_early, gap, intersect_func=intersect_func)
      data[[index]] <- simple_merge_prob(data[[index]], indi, nwindow[index])
    }

    nodes <- length(data[[index]])
    iteration <- 1
    if(verbose) cat("Start_CrossMerging_full\n")
    a <- start_end_block(data[[index]])
    a_old <- NULL
    while((length(a_old)==0 || (nrow(a_old)!=nrow(a)) || prod(a_old==a)==0) && ( length(a_old)==0 || (nrow(a_old)-nrow(a))>min_reduction_cross)){
      a_old <- a
      if(verbose) cat(paste("Iteration", iteration, ":", nodes, "nodes\n"))
      data[[index]] <- crossmerge(data[[index]], indi, nwindow[index], a, intersect_func=intersect_func)
      data[[index]] <- simple_merge_prob(data[[index]], indi, nwindow[index])
      a <- start_end_block(data[[index]])
      nodes <- length(data[[index]])
      iteration <- iteration + 1
    }


    nodes <- length(data[[index]])
    iteration <- 1
    if(verbose) cat("Start_IgnoreSmall\n")
    while((iteration==1 || (nrow(a_old)!=nrow(a)) || prod(a_old==a)==0) && ( iteration==1 || (nrow(a_old)- nrow(a))>min_reduction_neglet)){
      a_old <- a
      if(verbose) cat(paste("Iteration", iteration, ":", nodes, "nodes\n"))
      data[[index]] <- ignore_small_nodes(data[[index]], indi, nwindow[index], node_min, gap, intersect_func=intersect_func)
      data[[index]] <- simple_merge_prob(data[[index]], indi, nwindow[index], intersect_func=intersect_func)
      data[[index]] <- crossmerge(data[[index]], indi, nwindow[index], intersect_func=intersect_func)
      data[[index]] <- simple_merge_prob(data[[index]], indi, nwindow[index], intersect_func=intersect_func)
      a<- start_end_block(data[[index]])
      nodes <- length(data[[index]])
      iteration <- iteration + 1
    }

    dataset[[index]] <- block_dataset_construction(blockinfo[[index]], indi=indi, nwindow=nwindow[index])

    partial_blocklist[[index]] <- identify_blocks(data[[index]], indi, nwindow[index], min_share[index], edge_min=edge_min, subgroups=subgroups,
                                                  consider_nodes=consider_nodes, consider_edge=consider_edge, min_per_subgroup=min_per_subgroup,
                                                  subgroup_exception=subgroup_exception,
                                                  intersect_func=intersect_func, consider_multi=consider_multi, multi_min=multi_min,
                                                  double_share=double_share, node_min=node_min)
  }

  present_data <- matrix(1, nrow=indi, ncol = nwindow)
  if(deletion_count){
    for(index in 1:length(blockinfo[[1]])){
      for(index2 in 1:length(blockinfo[[1]][[index]][[6]])){
        present_data[blockinfo[[1]][[index]][[5]][[index2]], index] <- mean(blockinfo[[1]][[index]][[6]][[index2]]!=na_value)
      }
    }
  }

  blocklist <- list()
  for(index in 1:ncluster){
    if(length(partial_blocklist[[index]])>0){
      for(index2 in 1:length(partial_blocklist[[index]])){
        blocklist[[length(blocklist)+1]] <- partial_blocklist[[index]][[index2]]
        blocklist[[length(blocklist)]][[12]] <- index
      }
    }
  }



  # STRINGENZ ZU DATASET (ncol/nrow)

  current_iteration <- max_iteration
  if(length(target_coverage)>0){
    current_iteration <- 1
    blocklist_start <- blocklist
    min_majorblock_count <- numeric(max_iteration)
    coverage_results <- numeric(max_iteration)
    min_majorblock_count[1] <- min_majorblock
    blocklists <- list()

  }
  stop_iteration <- FALSE
  while(current_iteration <= max_iteration && !stop_iteration){

    if(length(target_coverage)>0){
      min_majorblock <- min_majorblock_count[current_iteration]
      blocklist <- blocklist_start
    }

    nodes <- length(blocklist)
    iteration <- 1
    if(verbose) cat("Start_Blockmerging\n")
    helper <- blocklist_startend(blocklist, type="snp")
    helper_old <- NULL
    while(iteration <= min_majorblock_steps || length(helper_old)==0 || (nrow(helper_old)!=nrow(helper)) || prod(helper_old==helper)==0){
      if(verbose) cat(paste("Iteration", iteration, ":", nodes, "blocks\n"))
      helper_old <- helper
      blocklist <- block_merging(blocklist, blockinfo, dataset, dhm, indi, nwindow, window_sequence_list, off_lines, min_similarity=min_similarity,
                                 consider_all=consider_all, node_min=node_min, save_allblock=save_allblock, helper=helper,
                                 c_dhm=c_dhm, c_dhm_mode=c_dhm_mode, intersect_func=intersect_func,
                                 min_per_subgroup=min_per_subgroup, subgroup_exception=subgroup_exception,
                                 run=(iteration-1))

      if(merge_closeblock==TRUE){
        blocklist <- block_closeblock_merging(blocklist, blockinfo, indi, nwindow, max_diff_l, max_diff_i, intersect_func=intersect_func,
                                              dataset=dataset)
      }
      if(min_majorblock>(0)){
        if(min_majorblock_steps>1){
          blocklist <- blockinfo_biggest(blocklist, nwindow, indi, min_majorblock=min(iteration-1, min_majorblock_steps-1)/(min_majorblock_steps-1)*min_majorblock,
                                         weighting_length=weighting_length, weighting_size=weighting_size,
                                         recalculate_biggest=recalculate_biggest, window_size=window_size,
                                         deletion_count=deletion_count, present_data=present_data)
        } else{
          blocklist <- blockinfo_biggest(blocklist, nwindow, indi, min_majorblock=min_majorblock, weighting_length=weighting_length, weighting_size=weighting_size,
                                         recalculate_biggest=recalculate_biggest, window_size=window_size,
                                         deletion_count=deletion_count, present_data=present_data)
        }
        if(length(blocklist)==0){

          if(length(target_coverage)>0){
            blocklist <- blocklist_start
            min_majorblock_count[current_iteration] <- ceiling(min(iteration-1, min_majorblock_steps-1)/(min_majorblock_steps-1)*min_majorblock / 2)
            min_majorblock <- min_majorblock_count[current_iteration]

          } else{
            if(verbose) cat(paste0("Empty Blocklist! Min_majorblock was chosen to high. Automatically set to ", ceiling(min(iteration-1, min_majorblock_steps-1)/(min_majorblock_steps-1)*min_majorblock / 5),"!\n"))
            if(verbose) cat("This might still be way too high for your data!!! \n")
            blocklist <- list()
            for(index in 1:ncluster){
              if(length(partial_blocklist[[index]])>0){
                for(index2 in 1:length(partial_blocklist[[index]])){
                  blocklist[[length(blocklist)+1]] <- partial_blocklist[[index]][[index2]]
                  blocklist[[length(blocklist)]][[12]] <- index
                }
              }
            }
            min_majorblock <- ceiling(min(iteration-1, min_majorblock_steps-1)/(min_majorblock_steps-1)*min_majorblock / 5)
          }
          if(min_majorblock<=1){
            stop("Empty blocklist with infinitisimal min_majorblock. Check input dataset.")
          }
          helper_old <- NULL
          iteration <- 0

        }
      }
      helper <- blocklist_startend(blocklist, type="snp")
      nodes <- length(blocklist)
      iteration <- iteration + 1
    }

    if(block_extending==TRUE){
      nodes <- length(blocklist)
      iteration <- 1
      extensions_done <- 0
      if(verbose) cat("Start_Blockextending\n")
      while(iteration==1 || (nrow(helper_old)!=nrow(helper)) || prod(helper_old==helper)==0 || extensions_done > 0){
        helper_old <- helper
        if(verbose) cat(paste("Iteration", iteration, ":", nodes, "blocks; ", extensions_done, "block extensions\n"))
        blocklist_out <- extend_block(blocklist, indi, nwindow, max_extending_diff=max_extending_diff,
                                      extending_ratio=extending_ratio, dataset=dataset, window_sequence_list=window_sequence_list)
        blocklist <- blocklist_out[[1]]
        extensions_done <- blocklist_out[[2]]


        blocklist <- block_merging(blocklist, blockinfo, dataset, dhm, indi, nwindow, window_sequence_list, off_lines, min_similarity=min_similarity,
                                   consider_all=consider_all, node_min=node_min, save_allblock=save_allblock, c_dhm=c_dhm,
                                   c_dhm_mode=c_dhm_mode , intersect_func=intersect_func,
                                   min_per_subgroup=min_per_subgroup, subgroup_exception=subgroup_exception,
                                   run=(iteration-1))

        if(length(blocklist)==0){

          if(length(target_coverage)>0){
            blocklist <- blocklist_start
            min_majorblock_count[current_iteration] <- ceiling(min(iteration-1, min_majorblock_steps-1)/(min_majorblock_steps-1)*min_majorblock / 2)
            min_majorblock <- min_majorblock_count[current_iteration]

          } else{
            if(verbose) cat(paste0("Empty Blocklist! Min_majorblock was chosen to high. Automatically set to ", ceiling(min(iteration-1, min_majorblock_steps-1)/(min_majorblock_steps-1)*min_majorblock / 5),"!\n"))
            if(verbose) cat("This might still be way too high for your data!!! \n")
            blocklist <- list()
            for(index in 1:ncluster){
              if(length(partial_blocklist[[index]])>0){
                for(index2 in 1:length(partial_blocklist[[index]])){
                  blocklist[[length(blocklist)+1]] <- partial_blocklist[[index]][[index2]]
                  blocklist[[length(blocklist)]][[12]] <- index
                }
              }
            }
            min_majorblock <- ceiling(min(iteration-1, min_majorblock_steps-1)/(min_majorblock_steps-1)*min_majorblock / 5)
          }
          if(min_majorblock<=1){
            stop("Empty blocklist with infinitisimal min_majorblock. Check input dataset.")
          }
          helper_old <- NULL
          iteration <- 0

        }


        if(merge_closeblock==TRUE){
          blocklist <- block_closeblock_merging(blocklist, blockinfo, indi, nwindow, max_diff_l, max_diff_i, intersect_func=intersect_func,
                                                dataset=dataset)
        }
        if(min_majorblock_steps>1){
          blocklist <- blockinfo_biggest(blocklist, nwindow, indi, min_majorblock=min(iteration-1, min_majorblock_steps-1)/(min_majorblock_steps-1)*min_majorblock,
                                         weighting_length=weighting_length, weighting_size=weighting_size,
                                         recalculate_biggest=recalculate_biggest, window_size=window_size,
                                         deletion_count=deletion_count, present_data=present_data)
        } else{
          blocklist <- blockinfo_biggest(blocklist, nwindow, indi, min_majorblock=min_majorblock,
                                         weighting_length=weighting_length, weighting_size=weighting_size,
                                         recalculate_biggest=recalculate_biggest, window_size=window_size,
                                         deletion_count=deletion_count, present_data=present_data)
        }

        helper <- blocklist_startend(blocklist, type="snp")
        nodes <- length(blocklist)
        iteration <- iteration + 1
      }
    }

    if(snp_extending==TRUE){
      # Fuege anfang und ende gleiche SNPs im Randbereich hinzu/ab
      blocklist <- extend_snp(blocklist, indi, nwindow, dhm, window_sequence_list, bp_map = bp_map,
                              max_extending_diff_snp=max_extending_diff_snp, extending_ratio_snp=extending_ratio_snp)

    }



    if(off_node_addition==TRUE){
      print("OFF_node_addition Only for single nwindow")
      blocklist <- add_offnode(blocklist=blocklist, dataset[[1]], indi, nwindow[[1]], window_sequence=window_sequence_list[[1]], bp_map=bp_map,
                               off_node_minimum_blocklength=off_node_minimum_blocklength,
                               off_node_minimum_blocksize=off_node_minimum_blocksize, raster=raster)
      nodes <- length(blocklist)
      iteration <- 1
      print("Start_Blockmerging")
      helper <- blocklist_startend(blocklist, type="snp")
      helper_old <- NULL
      while(length(helper_old)==0 || (nrow(helper_old)!=nrow(helper)) || prod(helper_old==helper)==0){
        print(paste("Iteration", iteration, ":", nodes, "blocks"))
        helper_old <- helper
        blocklist <- block_merging(blocklist, blockinfo, dataset, dhm, indi, nwindow, window_sequence_list, off_lines, min_similarity=min_similarity,
                                   consider_all=consider_all, node_min=node_min, save_allblock=save_allblock, helper=helper,
                                   c_dhm=c_dhm, c_dhm_mode=c_dhm_mode, intersect_func=intersect_func,
                                   min_per_subgroup=min_per_subgroup, subgroup_exception=subgroup_exception,
                                   run=(iteration-1))

        if(merge_closeblock==TRUE){
          blocklist <- block_closeblock_merging(blocklist, blockinfo, indi, nwindow, max_diff_l, max_diff_i, intersect_func=intersect_func,
                                                dataset=dataset)
        }
        if(min_majorblock>(0)){
          if(min_majorblock_steps>1){
            blocklist <- blockinfo_biggest(blocklist, nwindow, indi, min_majorblock=min(iteration-1, min_majorblock_steps-1)/(min_majorblock_steps-1)*min_majorblock,
                                           weighting_length=weighting_length, weighting_size=weighting_size,
                                           recalculate_biggest=recalculate_biggest, window_size=window_size,
                                           deletion_count=deletion_count, present_data=present_data)
          } else{
            blocklist <- blockinfo_biggest(blocklist, nwindow, indi, min_majorblock=min_majorblock, weighting_length=weighting_length, weighting_size=weighting_size,
                                           recalculate_biggest=recalculate_biggest, window_size=window_size,
                                           deletion_count=deletion_count, present_data=present_data)
          }
        }
        helper <- blocklist_startend(blocklist, type="snp")
        nodes <- length(blocklist)
        iteration <- iteration + 1
      }

      if(block_extending==TRUE){

        nodes <- length(blocklist)
        iteration <- 1
        extensions_done <- 0
        print("Start_Blockextending")
        while(iteration==1 || (nrow(helper_old)!=nrow(helper)) || prod(helper_old==helper)==0 || extensions_done > 0){
          helper_old <- helper
          print(paste("Iteration", iteration, ":", nodes, "blocks; ", extensions_done, "block extensions"))
          blocklist_out <- extend_block(blocklist, indi, nwindow, max_extending_diff=max_extending_diff,
                                        extending_ratio=extending_ratio, dataset=dataset, window_sequence_list=window_sequence_list)
          blocklist <- blocklist_out[[1]]
          extensions_done <- blocklist_out[[2]]

          blocklist <- block_merging(blocklist, blockinfo, dataset, dhm, indi, nwindow, window_sequence_list, off_lines, min_similarity=min_similarity,
                                     consider_all=consider_all, node_min=node_min, save_allblock=save_allblock, c_dhm=c_dhm,
                                     c_dhm_mode=c_dhm_mode , intersect_func=intersect_func,
                                     min_per_subgroup=min_per_subgroup, subgroup_exception=subgroup_exception,
                                     run=(iteration-1))

          if(merge_closeblock==TRUE){
            blocklist <- block_closeblock_merging(blocklist, blockinfo, indi, nwindow, max_diff_l, max_diff_i, intersect_func=intersect_func,
                                                  dataset=dataset)
          }
          if(min_majorblock_steps>1){
            blocklist <- blockinfo_biggest(blocklist, nwindow, indi, min_majorblock=min(iteration-1, min_majorblock_steps-1)/(min_majorblock_steps-1)*min_majorblock,
                                           weighting_length=weighting_length, weighting_size=weighting_size,
                                           recalculate_biggest=recalculate_biggest, window_size=window_size,
                                           deletion_count=deletion_count, present_data=present_data)
          } else{
            blocklist <- blockinfo_biggest(blocklist, nwindow, indi, min_majorblock=min_majorblock,
                                           weighting_length=weighting_length, weighting_size=weighting_size,
                                           recalculate_biggest=recalculate_biggest, window_size=window_size,
                                           deletion_count=deletion_count, present_data=present_data)
          }

          if(length(blocklist)==0){
            print("Value for min_majorblock to high! No blocks in haplotype library")
            return()
          }
          helper <- blocklist_startend(blocklist, type="snp")
          nodes <- length(blocklist)
          iteration <- iteration + 1
        }
      }

      if(snp_extending==TRUE){
        # Fuege anfang und ende gleiche SNPs im Randbereich hinzu/ab
        blocklist <- extend_snp(blocklist, indi, nwindow, dhm, window_sequence_list, bp_map = bp_map,
                                max_extending_diff_snp=max_extending_diff_snp, extending_ratio_snp=extending_ratio_snp)

      }


    }




    if(overlap_remove){

      if(verbose) cat("Start_Overlap_removal:\n")
      t <- coverage_test(blocklist, type="window")
      t1 <- coverage_test(blocklist, type="window", max=100)
      if(verbose) cat(paste0("Before: ", length(blocklist), " Blocks, ", round(100* mean(t), digits=2), " % Coverage, ", round(100 * (mean(t1)-mean(t)), digits=2)," % Overlapping segments.\n"))
      blocklist <- overlap_removal(blocklist, data, node_min=node_min)
      t <- coverage_test(blocklist, type="window")
      t1 <- coverage_test(blocklist, type="window", max=100)
      if(verbose) cat(paste0("After: ", length(blocklist), " Blocks, ", round(100* mean(t), digits=2), " % Coverage, ", round(100 * (mean(t1)-mean(t)), digits=2)," % Overlapping segments.\n"))

    }

    if(length(target_coverage)>0){

      blocklists[[current_iteration]] <- blocklist
      t <- coverage_test(blocklist, indi, type="snp")
      coverage_results[current_iteration] <- mean(t)
      prev_cov <- coverage_results[coverage_results>0]
      prev_block <- min_majorblock_count[min_majorblock_count>0]
      if(min(prev_cov) > target_coverage){
        min_majorblock_count[current_iteration + 1] <- max(max(prev_block)*2, max(prev_block)+min_step_size)
      } else if(max(prev_cov) < target_coverage){
        min_majorblock_count[current_iteration + 1] <- max(min(ceiling(min(prev_block)/2), min(prev_block) - min_step_size),1)
      } else{
        if(current_iteration==1){
          stop_iteration <- TRUE
        } else{
          pos <- ((prev_cov-target_coverage)>0) * (1:length(prev_cov))
          neg <- ((prev_cov-target_coverage)<0) * (1:length(prev_cov))
          if(FALSE){
            ordering <- sort(abs(prev_cov-target_coverage), index.return=TRUE)$ix
            neg <- base::intersect(ordering, neg)
            pos <- base::intersect(ordering, pos)
            min_majorblock_count[current_iteration + 1] <- ceiling(mean(prev_block[c(pos[1], neg[1])]))
          } else{
            min_majorblock_count[current_iteration +1 ] <- ceiling(mean(c(min(prev_block[neg]), max(prev_block[pos]))))
          }


        }

      }
      print(paste("Finish Target_Coverage Iteration", current_iteration))
      print(paste("Used min_majorblock:", min_majorblock_count[current_iteration]))
      print(paste("Achieved Coverage:", coverage_results[current_iteration]))
      if(current_iteration<max_iteration){
        print(paste("Start next iteration using min_majorblock:", min_majorblock_count[current_iteration+1]))
      }

      if(abs(coverage_results[current_iteration] - target_coverage)<target_stop){
        stop_iteration <- TRUE
      }




    }
    current_iteration <- current_iteration + 1

  }

  if(length(target_coverage)>0){
    print(paste("Final Iteration using min_majorblock",  min_majorblock_count[current_iteration-1]))
    print(paste("Achieved Coverage:", coverage_results[current_iteration-1]))
    take <- which.max(-abs(coverage_results-target_coverage))[1]
    blocklist <- blocklists[[take]]
  }


  if(major_snp_calculation==TRUE){
    blocklist <- major_snp_calculation(blocklist, dhm, recoding=recoding)
  }


  if(!developer_mode){
    if(big_output){
      long_blocklist <- blocklist
    }
    for(index in 1:length(blocklist)){
      blocklist[[index]][[12]] <- NULL
      blocklist[[index]][[11]] <- NULL
      blocklist[[index]][[10]] <- NULL
      blocklist[[index]][[9]] <- NULL
      blocklist[[index]][[8]] <- NULL
    }
  }



  # Leere Blöcke - wesentlich kleinere Blöcke nicht direkt entfernen sondern verkleinerte versionen behalten.
  if(big_output && developer_mode){
    return(list(blocklist, dataset, data, blockinfo, indi, nwindow))
  } else if(big_output){
    return(list(blocklist, dataset, data, blockinfo, indi, nwindow, long_blocklist, present_data))
  } else{
    return(blocklist)
  }

}
