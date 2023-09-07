#' read multiple bed with annotated data as per input parameter
#'
#' @param marker marker definition used to label folder and files eg MUTATIONS, LESIONS
#' @param figure figures like hypo/hyper to built the data path
#' @param sample_group name of the population used to build the data path
#'
#' @return list of pivot by column identified with column Label and by Sample

#'
read_multiple_bed <- function(sample_group, marker, figure)
{
  ssEnv <- get_session_info()

  f <- paste0(marker,"_", figure, sep="")
  souceFolder <- dir_check_and_create(ssEnv$result_folderData, c(as.character(sample_group),f))

  if(as.character(marker)=="DELTAS" | as.character(marker)=="DELTAQ"
    | as.character(marker)=="DELTAR" | as.character(marker)=="DELTARQ"
    | as.character(marker)=="BETA")
  {
    col_names <- c("CHR", "START", "END","VALUE","SAMPLEID")
  }
  else
  {
    col_names <- c("CHR", "START", "END", "SAMPLEID")
  }

  fileName <-file_path_build(souceFolder,c("MULTIPLE",as.character(marker),as.character(figure)),"fst")

  if(file.exists(fileName))
  {
    multiple_bed <- fst::read_fst(fileName, as.data.table = T)
    colnames(multiple_bed) <- col_names
    # just read the multiple bed
    if(!plyr::empty(multiple_bed))
    {
      multiple_bed[is.na(multiple_bed)] <- 0
      if(marker!="BETA" & marker!="DELTAR" & marker!="DELTAS" & marker!="DELTAQ" & !plyr::empty(multiple_bed))
        multiple_bed$VALUE <- 1

      multiple_bed <- data.frame(multiple_bed,"FIGURE" = as.character(figure), "MARKER" = as.character(marker), "SAMPLE_GROUP" = as.character(sample_group))
      multiple_bed$SAMPLE_GROUP <- as.factor(multiple_bed$SAMPLE_GROUP)
      multiple_bed$MARKER <- as.factor(multiple_bed$MARKER)
      multiple_bed$FIGURE <- as.factor(multiple_bed$FIGURE)
      return(multiple_bed)
    }
  }
}
