create_position_pivot_from_single_bed <- function(sample_group,sample_id,marker, figure, area, subarea){

  ssEnv <- get_session_info()

  marker <- as.character(marker)
  figure <- as.character(figure)
  area <- as.character(area)
  subarea <- as.character(subarea)

  pivot_file_name <- pivot_file_name(marker,figure,area,subarea)

  # read bed file
  bed_file <- bed_file_name(sample_id,sample_group,marker,figure)
  if (!file.exists(bed_file)){
    return()
  }

  data <- utils::read.delim(bed_file, sep="\t", header = FALSE, row.names = NULL, stringsAsFactors = FALSE)
  colnames(data) <- c("CHR","START","END",sample_id)

  # drop rows with NA in CHR
  data <- data[!is.na(data$START),]
  data <- data[!is.na(data$END),]
  data <- data[!is.na(data$CHR),]

  # look if pivot exists
  if (file.exists(pivot_file_name)){
    pivot <- readr::read_delim(pivot_file_name, show_col_types=FALSE, progress=FALSE,n_max=1)
    if(any(colnames(pivot) == sample_id))
      return()
    # read pivot as number except the first column
    pivot <- readr::read_delim(pivot_file_name,
      col_types = readr::cols(
        .default = readr::col_double(),
        CHR = readr::col_character(),
        START = readr::col_integer(),
        END = readr::col_integer()
        ),
      show_col_types=FALSE, progress=FALSE)
    if(marker=="MUTATIONS" || marker=="LESIONS")
      data[,sample_id] <- 1
    pivot <- merge(pivot, data, by = c("CHR", "START", "END"), all = TRUE)
  } else
  {
    pivot <- data
  }

  # save pivot
  readr::write_delim(pivot, pivot_file_name, col_names = T, progress = FALSE)
}
