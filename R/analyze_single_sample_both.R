analyze_single_sample_both <- function(envir, sample_detail) {

  start_time_single_sample <- Sys.time()
  # message(sample_detail$Sample_ID, " ", "SingleSampleBoth Sample analysis warmingUP ", Sys.time())
  result <- ""
  result <- result[-1]

  mutations <- data.frame("CHR"="", "START"="", "END"="")
  lesions <- mutations

  for( i in 1:nrow(envir$keys_figures))
  {
    figure <- envir$keys_figures[i,1]
    # if(figure=="BOTH")
    #   next
    folder_to_save <- dir_check_and_create(envir$result_folderData,c(as.character(sample_detail$Sample_Group),paste0("MUTATIONS","_", figure, sep = "")))
    fileName = file_path_build(folder_to_save,c(sample_detail$Sample_ID,"MUTATIONS",figure),"bed")
    if(file.exists(fileName))
    {
      # browser()
      mutationsTemp <- utils::read.csv(fileName, sep="\t", col.names =c("CHR", "START", "END") )
      mutations <- rbind(mutations, mutationsTemp )
    }

    folder_to_save <- dir_check_and_create(envir$result_folderData,c(as.character(sample_detail$Sample_Group),paste0("LESIONS","_", figure, sep = "")))
    fileName = file_path_build(folder_to_save,c(sample_detail$Sample_ID,"LESIONS",figure),"bed")
    if(file.exists(fileName))
    {
      # browser()
      lesionsTemp <- utils::read.csv(fileName, sep="\t", col.names =c("CHR", "START", "END") )
      lesions <- rbind(lesions, lesionsTemp )
    }

  }
  mutations <- mutations[-1,]
  lesions <- lesions[-1,]

  figure <- "BOTH"
  folder_to_save <- dir_check_and_create(envir$result_folderData,c(as.character(sample_detail$Sample_Group),paste0("MUTATIONS","_", figure, sep = "")))
  fileName = file_path_build(folder_to_save,c(sample_detail$Sample_ID,"MUTATIONS",figure),"bed")
  # mutations$SAMPLEID <- mutations$Sample_ID
  dump_sample_as_bed_file(
    data_to_dump = mutations,
    fileName = fileName
  )

  # browser()
  folder_to_save <- dir_check_and_create(envir$result_folderData,c(as.character(sample_detail$Sample_Group),paste0("LESIONS","_", figure, sep = "")))
  fileName = file_path_build(folder_to_save,c(sample_detail$Sample_ID,"LESIONS",figure),"bed")
  lesions <- unique(lesions)
  if(nrow(lesions)>0)
  {
    # lesions$SAMPLEID <- sample_detail$Sample_ID
    dump_sample_as_bed_file(
      data_to_dump = lesions,
      fileName = fileName
    )
  }

  result["LESIONS_BOTH"] <- if (!is.null(lesions)) nrow(lesions) else 0
  result["MUTATIONS_BOTH"] <- if (!is.null(mutations)) nrow(mutations) else 0

  end_time_single_sample <- Sys.time()
  time_taken <- end_time_single_sample - start_time_single_sample
  # message(sample_detail$Sample_ID, " ", "Completed sample ", time_taken)

  return(result)
}
