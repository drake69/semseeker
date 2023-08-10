read_annotated_bed <- function(
  figure,
  marker,
  area,
  subarea)
{

  figure <- as.character(figure)
  marker <- as.character(marker)
  area <- as.character(area)
  subarea <- as.character(subarea)

  ssEnv <- get_session_info()
  # area and subarea are defined using the filename
  dest_folder <- dir_check_and_create(ssEnv$result_folderData,subFolders = c("Annotated"))
  bedFileName <- file_path_build(dest_folder , c(marker, figure, area,subarea, "Annotated"),"fst")
  if(file.exists(bedFileName))
  {
    if(file.info(bedFileName)$size < 10)
    {
      message("WARNING: ", Sys.time(), " Given up file:", bedFileName, " is empty!")
    }
    else
    {
      final_bed <- fst::fst(path = bedFileName)
      final_bed$VALUE = as.numeric(final_bed$VALUE)
    }
    final_bed_temp <- as.data.frame(final_bed)
    return (final_bed_temp)
  }

  return(NULL)
}
