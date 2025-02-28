#' Title
#'
#' @param sample_sheet
#'
#' @returns summary of sample sheet
#'
#' @importFrom doRNG %dorng%
#'
deltaX_get <- function(sample_sheet)
{
  ssEnv <- get_session_info()
  # ssEnv <- get_session_info("~/Documents/Dati_Lavoro/cancer_stage/results/ewas_data_hub_B_Q2/")

  keys <- ssEnv$keys_markers_figures
  # remove figure colunn
  keys <- keys[,-which(colnames(keys)=="FIGURE")]
  keys <- keys[,-which(colnames(keys)=="COMBINED")]
  keys <- unique(keys)

  area_position <-"POSITION"
  subarea_position <- ""

  if(ssEnv$showprogress)
    progress_bar <- progressr::progressor(along = 1:nrow(keys))
  else
    progress_bar <- ""

  variables_to_export <- c("ssEnv", "dir_check_and_create", "subarea",
    "progress_bar","progression_index", "progression", "progressor_uuid",
    "owner_session_uuid", "trace","probe_features_get", "localKeys",
    "file_path_build","%>%","get_session_info","log_event")
  sample_sheet_res <- data.frame()
  keys <- subset(keys,!is.na(SOURCE))
  keys <- subset(keys,!is.na(Q))
  keys <- subset(keys,(Q!=1))
  if(nrow(keys)==0)
    return(sample_sheet)

  # sample_sheet_res <-  foreach::foreach(k =1:nrow(keys), .combine= cbind , .export = variables_to_export) %dorng%
  for ( k in 1:nrow(keys))
  {
    # browser()
    sample_sheet_temp <- as.data.frame(sample_sheet[,c("Sample_ID","Sample_Group")])

    key <- keys[k,]
    source_marker <- key$SOURCE
    figure <- key$FIGURE

    pivot_file_nameparquet <- pivot_file_name_parquet(source_marker,"HYPER",area_position,subarea_position)
    pivot_hyper <- polars::pl$scan_parquet(pivot_file_nameparquet)
    positions_hyper <- pivot_hyper$select(pivot_hyper$columns[1:3])$collect()$to_data_frame()
    pivot_hyper <- pivot_hyper$drop(pivot_hyper$columns[1:3])
    vector_shaped_hyper <- as.vector(as.matrix(pivot_hyper$collect()$to_data_frame()))
    dim_pivot_hyper <- dim(pivot_hyper$collect()$to_data_frame())
    colname_pivot_hyper <- colnames(pivot_hyper)
    rm(pivot_hyper)
    vector_shaped_hyper[vector_shaped_hyper==0] <- NA
    length_hyper <- length(vector_shaped_hyper)

    pivot_file_nameparquet <- pivot_file_name_parquet(source_marker,"HYPO",area_position,subarea_position)
    pivot_hypo <- polars::pl$scan_parquet(pivot_file_nameparquet)
    positions_hypo <- pivot_hypo$select(pivot_hypo$columns[1:3])$collect()$to_data_frame()
    pivot_hypo <- pivot_hypo$drop(pivot_hypo$columns[1:3])
    vector_shaped_hypo <- as.vector(as.matrix(pivot_hypo$collect()$to_data_frame()))
    dim_pivot_hypo <- dim(pivot_hypo$collect()$to_data_frame())
    colname_pivot_hypo <- colnames(pivot_hypo)
    rm(pivot_hypo)
    vector_shaped_hypo[vector_shaped_hypo==0] <- NA
    length_hypo <- length(vector_shaped_hypo)

    vector_shaped <- c(vector_shaped_hyper,vector_shaped_hypo)
    positions <- rbind(positions_hyper,positions_hypo)

    ####
    marker <- key$MARKER
    if(endsWith(marker,"P"))
      vector_both <-  cut(vector_shaped, breaks=as.numeric(key$Q), labels=FALSE)
    else if(endsWith(marker,"Q"))
      vector_both <- as.numeric(dplyr::ntile(x=vector_shaped , n= as.numeric(key$Q)))

    sample_sheet_temp <- save_figure(colname_pivot_hyper, dim_pivot_hyper,vector_both[(1:length_hyper)],positions_hyper,sample_sheet_temp,area_position,subarea_position,marker,"HYPER")
    sample_sheet_temp <- save_figure(colname_pivot_hypo, dim_pivot_hypo,vector_both[((length_hyper + 1):(length_hyper + length_hypo))], positions_hypo,sample_sheet_temp,area_position,subarea_position,marker,"HYPO")

    # if ( nrow(sample_sheet_res)==0)
    #   sample_sheet_res <- sample_sheet_temp
    # else
    #   sample_sheet_res <- merge(sample_sheet_res,sample_sheet_temp,by="Sample_ID")

    if(ssEnv$showprogress)
      progress_bar(sprintf("Got: %s", stringr::str_pad(marker, 10, pad = " ")))

    # remove Sample_ID column
    sample_sheet_temp <- sample_sheet_temp[,-which(colnames(sample_sheet_temp)=="Sample_ID")]
    sample_sheet_temp <- sample_sheet_temp[,-which(colnames(sample_sheet_temp)=="Sample_Group")]
    sample_sheet_temp
  }
  if ( nrow(sample_sheet_res)==nrow(sample_sheet))
    sample_sheet <- cbind(sample_sheet,sample_sheet_res)
  return(sample_sheet)
}


save_figure <- function(colname_pivot,dim_pivot,vector_figure,positions,sample_sheet,area, subarea,marker,figure)
{
  vector_figure <- matrix(vector_figure, nrow=dim_pivot[1], ncol=dim_pivot[2])
  vector_figure <- as.data.frame(vector_figure)
  colnames(vector_figure) <- colname_pivot
  vector_figure <- cbind(positions,vector_figure)
  pivot_file_nameparquet <- pivot_file_name_parquet(marker,figure,area,subarea)
  vector_figure <- polars::as_polars_df(as.data.frame(vector_figure))
  vector_figure <- vector_figure$sort(c("CHR", "START"), descending = c(FALSE, FALSE))
  vector_figure$write_parquet(pivot_file_nameparquet)

  vector_figure <- vector_figure$to_data_frame()
  marker_sum <- data.frame(colSums(vector_figure[,4:ncol(vector_figure)],na.rm = T))
  colnames(marker_sum) <- paste0(marker,"_",figure)
  marker_sum$Sample_ID <- row.names(marker_sum)
  sample_sheet <- merge(sample_sheet,marker_sum,by="Sample_ID",all.x=TRUE)

  if(ssEnv$showprogress)
    progress_bar <- progressr::progressor(along = 4:ncol(vector_figure))
  else
    progress_bar <- ""

  for(c in 4:(ncol(vector_figure)))
  {
    vector_figure_sample <- vector_figure[,c]
    sample_id <- colnames(vector_figure)[c]
    sample_group <- sample_sheet[sample_sheet$Sample_ID==sample_id,"Sample_Group"]
    result <- cbind(positions,vector_figure_sample)
    # drop rows with NA
    result <- result[complete.cases(result),]
    result <- result[result[,4]>0,]
    bed_file_name <- bed_file_name(sample_id,sample_group,marker,figure)
    dump_sample_as_bed_file(result,bed_file_name)
    if(ssEnv$showprogress)
      progress_bar(sprintf("Annotating bed file %s of %s, %s.",sample_id,marker, figure))
  }
  # replace NA with 0
  sample_sheet[is.na(sample_sheet)] <- 0
  return(sample_sheet)
}
