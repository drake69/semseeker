get_pivot_both <- function(marker)
{
  # marker <- "DELTAP"
  # ssEnv <- get_session_info("~/Documents/Dati_Lavoro/osteoporosis/results/GSE99624")
  # update_session_info(ssEnv)
  ssEnv <- get_session_info()
  is_discrete <- unique(ssEnv$keys_markers_figures_default[ssEnv$keys_markers_figures_default$MARKER==marker,"DISCRETE"])

  pivot_file_name_both <- pivot_file_name_parquet(marker,"BOTH","PROBE","")
  if(file.exists(pivot_file_name_both))
    return(arrow::read_parquet(pivot_file_name_both))

  pivot_file_name_hyper <- pivot_file_name_parquet(marker,"HYPER","PROBE","")
  pivot_file_name_hypo <- pivot_file_name_parquet(marker,"HYPO","PROBE","")

  if(!file.exists(pivot_file_name_hyper) && !file.exists(pivot_file_name_hypo))
    return(data.frame())
  pivot_hyper <- polars::pl$scan_parquet(pivot_file_name_hyper)
  pivot_hypo <- polars::pl$scan_parquet(pivot_file_name_hypo)

  # Union (concatenate) the two DataFrames
  pivot_both <- polars::pl$concat(list(pivot_hyper, pivot_hypo))$collect()
  pivot_both <- pivot_both$group_by("AREA", maintain_order=FALSE)$sum()

  pivot_both$write_parquet(pivot_file_name_both)

  # pivot_both <- pivot_both$group_by("AREA", maintain_order=FALSE)$agg(
  #   pl$all()$exclude("AREA")$sum()
  # )
  #
  # pivot_both$sink_parquet(pivot_file_name_both)


  # pivot_hyper <- data.frame()
  # if(file.exists(pivot_file_name_hyper))
  #   pivot_hyper <- polars::pl$read_parquet(pivot_file_name_hyper)$to_data_frame()
  # pivot_hypo <- data.frame()
  # if(file.exists(pivot_file_name_hypo))
  #   pivot_hypo <- polars::pl$read_parquet(pivot_file_name_hypo)$to_data_frame()
  #
  # # get the row with the max count of values
  # count_m <- plyr::rbind.fill(as.data.frame(pivot_hyper), as.data.frame(pivot_hypo))
  # rm(pivot_hypo)
  # rm(pivot_hyper)
  #
  #
  # # sort count_m by AREA
  # count_m <- count_m[order(count_m$AREA),]
  #
  # if (nrow(count_m)!=0)
  #   # Group by AREA and sum all other columns
  #   df_grouped <- df$group_by("AREA")$agg(
  #     pl$all()$exclude("AREA")$sum()
  #   )
  #
  # count_m <- as.data.frame(count_m)
  # arrow::write_parquet(count_m, pivot_file_name_both)
  return(pivot_both$to_data_frame())
}
