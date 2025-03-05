#' calculate scores for a study across all the used quantile/bin number with result saved into a base folder
#' @importFrom doRNG %dorng%
#' @importFrom doFuture %dofuture%
#' @export
marker_quantization_metric_q_search <- function(study, qs=c(2,4,8,16,32), base_folder,dest_folder,folder_key,
  maxResources =  maxResources, parallel_strategy  =  parallel_strategy, start_fresh = FALSE, ...)
{

  to_export <- ""
  result_folder <- paste0(base_folder,"/", study, folder_key, qs[1], sep="")
  ssEnv <- init_env( result_folder =  result_folder, maxResources =  maxResources, parallel_strategy  =  parallel_strategy,
    start_fresh = FALSE, areas=c("POSITION"), subareas=c(""), markers=c("DELTAS","DELTAR"), figures=c("HYPO","HYPER"), ...)

  study_summary <-   study_summary_get()
  create_position_pivots(study_summary,ssEnv$keys_areas_subareas_markers_figures)

  nkeys <- 2*2*length(qs) + 1
  if(ssEnv$showprogress)
    progress_bar <- progressr::progressor(along = 1:(nkeys))
  else
    progress_bar <- ""
  log_event("BANNER:", format(Sys.time(), "%a %b %d %X %Y"), " SemSeeker will evaluate the best quantile/bin as deltas/deltar discretisation for project \n in ", ssEnv$result_folderData)


  dataFolder <- dir_check_and_create(dest_folder,c("Data/Distributions"))
  filename  =  file_path_build(dataFolder,c(study,"Q_SEARCH_DISTRIBUTION_ANALYSIS_ALL"),"csv")
  # message(filename)
  result_temp <- data.frame()
  if(file.exists(filename))
    result_temp <- utils::read.csv2(file = filename)

  qq <- nrow(result_temp)==0
  if(!qq)
      qq <- nrow(subset(result_temp,MARKER=="MUTATIONS"))==0

  if(qq)
  {
    res_temp <- data.frame("ORIGINAL_MARKER"="DELTAS")
    original <- load_deltax("DELTAS")
    log_event("DEBUG:", format(Sys.time(), "%a %b %d %X %Y"), "Processing key: MUTATIONS")
    res_temp$MARKER <- "MUTATIONS"
    quantized <- rep(1, length(original))
    mdl_perf <- model_performance(original, quantized,c(),c())
    res_temp <- cbind(res_temp, mdl_perf)
    # SSIM con pracma
    ssim_value <- ssim(original, quantized)
    res_temp$Structural_Similarity_Index <- ssim_value
    vi_value <- variation_of_information(original, quantized)
    res_temp$variation_of_information <- vi_value
    jsd <- jsd_calc(original, quantized)
    res_temp$JSD <- jsd
    result_temp <- data.frame(res_temp,"Q"=qs)
  }
  if(ssEnv$showprogress)
    progress_bar(sprintf("Doing comparison."))

  for (q in qs)
  {
    for ( k in c("DELTAS","DELTAR"))
    {
      original <- load_deltax(k)
      for(ss in c("Q","P"))
      {
        res_temp <- data.frame("Q"=q)
        res_temp$ORIGINAL_MARKER <- k

        if(ss=="Q" & k=="DELTAS")
          res_temp$MARKER <- "DELTAQ"

        if(ss=="Q" & k=="DELTAR")
          res_temp$MARKER <- "DELTARQ"

        if(ss=="P" & k=="DELTAS")
          res_temp$MARKER <- "DELTAP"

        if(ss=="P" & k=="DELTAR")
          res_temp$MARKER <- "DELTARP"

        # skip already processed
        if(nrow(result_temp)>0)
          it_exists <- nrow(subset(result_temp,(MARKER == res_temp$MARKER) & (Q == q))) > 0
        else
          it_exists <- FALSE
        if(it_exists)
        {
          if(ssEnv$showprogress)
            progress_bar(sprintf("Doing comparison."))
          next
        }

        log_event("DEBUG:", format(Sys.time(), "%a %b %d %X %Y"), "Processing key: ", res_temp$MARKER, "for quantile:", q)

        if(ss=="Q")
          quantized <- as.numeric(dplyr::ntile(x=original , n=as.numeric(q)))

        if(ss=="P")
          quantized <- cut(original, breaks=q, labels=FALSE)

        quantized <- quantized * (max(original) / q)

        mdl_perf <- model_performance(original, quantized,c(),c())
        res_temp <- cbind(res_temp, mdl_perf)

        # SSIM con pracma
        ssim_value <- ssim(original, quantized)
        res_temp$Structural_Similarity_Index <- ssim_value

        vi_value <- variation_of_information(original, quantized)
        res_temp$variation_of_information <- vi_value

        jsd <- jsd_calc(original, quantized)
        res_temp$JSD <- jsd

        result_temp <- plyr::rbind.fill(result_temp, res_temp)
        if(ssEnv$showprogress)
          progress_bar(sprintf("Doing comparison."))

      }
      # res_temp
    }
  }



  #
  if(nrow(result_temp)==0)
  {
    log_event("DEBUG:", format(Sys.time(), "%a %b %d %X %Y")," Empty result_temp")
    return(NULL)
  }


  ssEnv <- init_env( result_folder =  dest_folder, maxResources =  maxResources, parallel_strategy  =  parallel_strategy,
    start_fresh = FALSE, areas=c("PROBE"), subareas=c("WHOLE"), ...)

  colnames(result_temp) <- toupper(colnames(result_temp))
  dataFolder <- dir_check_and_create(dest_folder,c("Data","Distributions"))
  filename  <-  file_path_build(dataFolder,c(study,"Q_SEARCH_DISTRIBUTION_ANALYSIS_ALL"),"csv")
  utils::write.csv2(result_temp, file = filename, row.names = FALSE)

  markers <- unique(result_temp$MARKER)
  filename  <-  file_path_build(dest_folder,c(study,"Q_SEARCH_DISTRIBUTION_ANALYSIS_PIVOT"),"csv")
  aggregate_pivot <- data.frame()
  for (m in 1:length(markers))
  {

    # m <- 1
    marker <- markers[m]
    result_temp_marker <- result_temp[result_temp$MARKER == marker,]

    # cols_to_cycke <- colnames(result_temp_marker)[4:15]
    cols_to_cycke <- ssEnv$model_metrics
    cols_to_cycke <- which(colnames(result_temp_marker) %in% cols_to_cycke)
    scores <- data.frame()
    for(metric in cols_to_cycke)
    {
      metric <- colnames(result_temp_marker)[metric]
      scores_temp <- data.frame()
      scores_temp <- metrics_ranking(metric,result_temp_marker,column_to_rank =metric)
      scores_temp$FIGURE <- result_temp_marker$FIGURE
      scores <- plyr::rbind.fill(scores, scores_temp)
    }

    scores$SCORE <- round(scores$SCORE, 2)
    # aggregate scores by MARKER and sum RANK
    # scores_agg <- aggregate(scores$SCORE, by = list(scores$Q), FUN = sum)
    # calculate TOTAL per marker
    # scores_agg <- scores_agg[order(scores_agg$x, decreasing = TRUE),]
    scores <- scores[,c("Q","SCORE")]

    # # create a pivot table with Q as rows and FIGURE as columns
    # scores_agg_fig <- merge(scores_agg_fig,scores_agg, by="Q")
    # # sort by SCORE descending
    # scores_agg_fig <- scores_agg_fig[order(scores_agg_fig$TOTAL, decreasing = TRUE),]
    scores$MARKER <- marker
    # save scores
    aggregate_pivot <- plyr::rbind.fill(aggregate_pivot, scores)
  }
  library(dplyr)
  filename  =  file_path_build(dataFolder,c(study,"Q_SEARCH_DISTRIBUTION", "ANALYSIS","SCORE"),"csv")
  # sum the score over Q and marker
  result <- aggregate_pivot %>%
    group_by(MARKER, Q) %>%
    dplyr::summarise(SCORE = sum(SCORE, na.rm = TRUE))

  utils::write.csv2(result, file = filename, row.names = FALSE)
  aggregate_pivot <- reshape2::dcast(aggregate_pivot, MARKER ~ Q, value.var = "SCORE", fun.aggregate = sum)
  aggregate_pivot <- t(aggregate_pivot)
  aggregate_pivot <- as.data.frame(aggregate_pivot)
  colnames(aggregate_pivot) <- aggregate_pivot[1,]
  aggregate_pivot <- aggregate_pivot[-1,]
  aggregate_pivot$Q <- rownames(aggregate_pivot)
  filename  =  file_path_build(dataFolder,c(study,"Q_SEARCH_DISTRIBUTION", "ANALYSIS","SCORE","PIVOT"),"csv")
  utils::write.csv2(aggregate_pivot, file = filename, row.names = FALSE)


  # calculate best marker

  filename  =  file_path_build(dataFolder,c(study,"Q_SEARCH_DISTRIBUTION_ANALYSIS_ALL"),"csv")
  result_temp <- utils::read.csv2(file = filename)
  filename  =  file_path_build(dataFolder,c(study,"Q_SEARCH_DISTRIBUTION", "ANALYSIS","SCORE"),"csv")
  results_q <- utils::read.csv2(file =filename)
  results_q_max <- unique(results_q[,c("SCORE","MARKER")] %>%
      dplyr::group_by(MARKER) %>%
      dplyr::filter(SCORE == max(SCORE)))

  # get Q corresponding to max score
  results_q_max <- merge(results_q, results_q_max, by = c("MARKER","SCORE"))

  # remove mutations
  results_q_max <- results_q_max[results_q_max$MARKER != "MUTATIONS",]

  # add MUTATIONS
  results_q_max <- rbind(results_q_max, data.frame("MARKER"="MUTATIONS","SCORE"=5,"Q"=2))

  # filter results_q with results_q_max by MARKER and SCORE
  results_q <- merge(result_temp, results_q_max, by = c("MARKER","Q"))

  # anaysis cross marker
  cols_to_cycke <- ssEnv$model_metrics
  cols_to_cycke <- which(colnames(results_q) %in% cols_to_cycke)
  scores <- data.frame()
  for(metric in cols_to_cycke)
  {
    metric <- colnames(results_q)[metric]
    scores_temp <- data.frame()
    scores_temp <- metrics_ranking(metric,results_q,column_to_rank =metric)
    scores <- plyr::rbind.fill(scores, scores_temp)
  }

  scores$SCORE <- round(scores$SCORE, 2)
  scores <- scores[,c("MARKER","SCORE")]
  # sum scores per MARKER
  scores <- scores %>%
    dplyr::group_by(MARKER) %>%
    dplyr::summarise(SCORE = sum(SCORE))
  filename  =  file_path_build(dataFolder,c(study,"BEST_MARKER", "ANALYSIS","SCORE"),"csv")
  # sort by score desc
  scores <- scores[order(-scores$SCORE),]
  utils::write.csv2(scores, file = filename, row.names = FALSE)

  close_env()
}


load_deltax <- function(source_marker){

  ssEnv <- get_session_info()
  area_position <- "POSITION"
  subarea_position <- ""
  log_event("DEBUG: ", format(Sys.time(), "%a %b %d %X %Y"), " loading: ", source_marker)

  pivot_file_nameparquet <- pivot_file_name_parquet(source_marker,"HYPER",area_position,subarea_position)
  pivot_hyper <- polars::pl$scan_parquet(pivot_file_nameparquet)
  positions_hyper <- pivot_hyper$select(pivot_hyper$columns[1:3])$collect()$to_data_frame()
  pivot_hyper <- pivot_hyper$drop(pivot_hyper$columns[1:3])
  vector_shaped_hyper <- as.vector(as.matrix(pivot_hyper$collect()$to_data_frame()))
  rm(pivot_hyper)
  vector_shaped_hyper[vector_shaped_hyper==0] <- NA

  pivot_file_nameparquet <- pivot_file_name_parquet(source_marker,"HYPO",area_position,subarea_position)
  pivot_hypo <- polars::pl$scan_parquet(pivot_file_nameparquet)
  positions_hypo <- pivot_hypo$select(pivot_hypo$columns[1:3])$collect()$to_data_frame()
  pivot_hypo <- pivot_hypo$drop(pivot_hypo$columns[1:3])
  vector_shaped_hypo <- as.vector(as.matrix(pivot_hypo$collect()$to_data_frame()))
  rm(pivot_hypo)
  vector_shaped_hypo[vector_shaped_hypo==0] <- NA

  vector_shaped <- c(vector_shaped_hyper,vector_shaped_hypo)
  vector_shaped <- vector_shaped[!is.na(vector_shaped)]
  log_event("DEBUG: ", format(Sys.time(), "%a %b %d %X %Y"), " loaded ", source_marker)
  return(vector_shaped)
}

jsd_calc <- function(original, quantized){
  # calculate JSD
  common_events <- unique(c(original, quantized))
  frequency_table1_adjusted <- tabulate(match(original, common_events), nbins = length(common_events))
  frequency_table2_adjusted <- tabulate(match(quantized, common_events), nbins = length(common_events))
  probability_distribution1_adjusted <- frequency_table1_adjusted / sum(frequency_table1_adjusted)
  probability_distribution2_adjusted <- frequency_table2_adjusted / sum(frequency_table2_adjusted)
  jsd <- suppressMessages(suppressWarnings(philentropy::JSD(rbind(probability_distribution1_adjusted, probability_distribution2_adjusted))))
  return(jsd)

}
