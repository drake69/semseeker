deltaX_get <- function(sample_sheet)
{
  ssEnv <- get_session_info()
  keys <- ssEnv$keys_markers_figures
  # remove figure colunn
  keys <- keys[,-which(colnames(keys)=="FIGURE")]
  keys <- keys[,-which(colnames(keys)=="COMBINED")]
  keys <- unique(keys)

  area_position <-"POSITION"
  subarea_position <- ""

  for ( k in 1:nrow(keys))
  {
    key <- keys[k,]
    if(is.na(key$SOURCE))
      next
    # if(is.na(key$Q))
    if(key$Q==1 | is.na(key$Q))
      next
    source_marker <- key$SOURCE
    figure <- key$FIGURE
    pivot_file_name <- pivot_file_name(source_marker,"HYPER",area_position,subarea_position)
    pivot_hyper <- readr::read_delim(pivot_file_name,
      col_types = vroom::cols(
        .default = vroom::col_double(),
        CHR = vroom::col_character(),
        START = vroom::col_integer(),
        END = vroom::col_integer()
      ),
      show_col_types=FALSE, progress=FALSE)
    positions_hyper <- pivot_hyper[,c(1,2,3)]
    pivot_hyper <- pivot_hyper[,-c(1,2,3)]
    vector_shaped_hyper <- as.vector(as.matrix(pivot_hyper))
    vector_shaped_hyper[vector_shaped_hyper==0] <- NA

    pivot_file_name <- pivot_file_name(source_marker,"HYPO",area_position,subarea_position)
    pivot_hypo <- readr::read_delim(pivot_file_name,
      col_types = vroom::cols(
        .default = vroom::col_double(),
        CHR = vroom::col_character(),
        START = vroom::col_integer(),
        END = vroom::col_integer()
      ),
      show_col_types=FALSE, progress=FALSE)
    # align colnames

    positions_hypo <- pivot_hypo[,c(1,2,3)]
    pivot_hypo <- pivot_hypo[,-c(1,2,3)]
    pivot_hypo <- pivot_hypo[,colnames(pivot_hyper)]
    vector_shaped_hypo <- as.vector(as.matrix(pivot_hypo))
    vector_shaped_hypo[vector_shaped_hypo==0] <- NA

    vector_shaped <- c(vector_shaped_hyper,vector_shaped_hypo)
    positions <- rbind(positions_hyper,positions_hypo)

    #### BOTH
    marker <- key$MARKER
    if(endsWith(marker,"P"))
      delta_x <-  cut(vector_shaped, breaks=as.numeric(key$Q), labels=FALSE)
    else if(endsWith(marker,"Q"))
      delta_x <- as.numeric(dplyr::ntile(x=vector_shaped , n= as.numeric(key$Q)))
    # else
    #   delta_x <- as.numeric(vector_shaped>0)

    delta_x[is.na(delta_x)] <- 0
    delta_x <- matrix(delta_x, nrow=nrow(pivot_hyper) + nrow(pivot_hypo), ncol=ncol(pivot_hyper))
    delta_x <- as.data.frame(delta_x)
    colnames(delta_x) <- colnames(pivot_hyper)

    result <- cbind(positions,delta_x)
    pivot_file_name <- pivot_file_name(marker,"BOTH",area_position,subarea_position)
    readr::write_delim(result,pivot_file_name, col_names = T, progress = FALSE, delim=";")

    marker_both <- data.frame(colSums(delta_x))
    colnames(marker_both) <- paste0(marker,"_","BOTH")
    marker_both$Sample_ID <- row.names(marker_both)
    sample_sheet <- merge(sample_sheet,marker_both,by.x="Sample_ID",all=TRUE)

    #### HYPER
    delta_x_hyper <- delta_x[1:nrow(pivot_hyper),]
    result <- cbind(positions_hyper,delta_x_hyper)
    pivot_file_name <- pivot_file_name(marker,"HYPER",area_position,subarea_position)
    readr::write_delim(result,pivot_file_name, col_names = T, progress = FALSE, delim=";")

    marker_hyper <- data.frame(colSums(delta_x_hyper))
    colnames(marker_hyper) <- paste0(marker,"_","HYPER")
    marker_hyper$Sample_ID <- row.names(marker_hyper)
    sample_sheet <- merge(sample_sheet,marker_hyper,by.x="Sample_ID",all=TRUE)

    for(c in 1:ncol(pivot_hyper))
    {
      delta_x_hyper_sample <- delta_x_hyper[,c]
      sample_id <- colnames(pivot_hyper)[c]
      sample_group <- sample_sheet[sample_sheet$Sample_ID==sample_id,"Sample_Group"]
      result <- cbind(positions_hyper,delta_x_hyper_sample)
      result <- result[result[,4]>0,]
      bed_file_name <- bed_file_name(sample_id,sample_group,marker,"HYPER")
      dump_sample_as_bed_file(result,bed_file_name)
    }

    #### HYPO
    delta_x_hypo <- delta_x[(nrow(pivot_hyper) + 1):nrow(delta_x),]
    result <- cbind(positions_hypo,delta_x_hypo)
    pivot_file_name <- pivot_file_name(marker,"HYPO",area_position,subarea_position)
    readr::write_delim(result,pivot_file_name, col_names = T, progress = FALSE, delim=";")

    marker_hypo <- data.frame(colSums(delta_x_hypo))
    colnames(marker_hypo) <- paste0(marker,"_","HYPO")
    marker_hypo$Sample_ID <- row.names(marker_hypo)
    sample_sheet <- merge(sample_sheet,marker_hypo,by.x="Sample_ID",all=TRUE)

    for(c in 1:ncol(pivot_hypo))
    {
      delta_x_hypo_sample <- delta_x_hypo[,c]
      sample_id <- colnames(pivot_hypo)[c]
      sample_group <- sample_sheet[sample_sheet$Sample_ID==sample_id,"Sample_Group"]
      result <- cbind(positions_hypo,delta_x_hypo_sample)
      result <- result[result[,4]>0,]
      bed_file_name <- bed_file_name(sample_id,sample_group,marker,"HYPO")
      dump_sample_as_bed_file(result,bed_file_name)
    }
  }
  return(sample_sheet)
}
