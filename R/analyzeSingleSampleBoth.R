analyzeSingleSampleBoth <- function(sampleDetail) {

  start_time_single_sample <- Sys.time()
  message(sampleDetail$Sample_ID, " ", "... Sample analysis warmingUP ", Sys.time())
  result <- ""
  result <- result[-1]

  mutations <- data.frame("CHR"="", "START"="", "END"="")
  lesions <- mutations

  for( figure in keys.figures)
  {
    if(figure=="BOTH")
      next
    folder2Save <- dir_check_and_create(resultFolderData,c(as.character(sampleDetail$Sample_Group),paste0("MUTATIONS","_", figure, sep = "")))
    fileName = file_path_build(folder2Save,c(sampleDetail$Sample_ID,"MUTATIONS",figure),"bed")
    if(file.exists(fileName))
    {
      # browser()
      mutationsTemp <- read.csv(fileName, sep="\t", col.names =c("CHR", "START", "END") )
      mutations <- rbind(mutations, mutationsTemp )
    }

    folder2Save <- dir_check_and_create(resultFolderData,c(as.character(sampleDetail$Sample_Group),paste0("LESIONS","_", figure, sep = "")))
    fileName = file_path_build(folder2Save,c(sampleDetail$Sample_ID,"LESIONS",figure),"bed")
    if(file.exists(fileName))
    {
      # browser()
      lesionsTemp <- read.csv(fileName, sep="\t", col.names =c("CHR", "START", "END") )
      lesions <- rbind(lesions, lesionsTemp )
    }

  }
  mutations <- mutations[-1,]
  lesions <- lesions[-1,]

  figure <- "BOTH"
  folder2Save <- dir_check_and_create(resultFolderData,c(as.character(sampleDetail$Sample_Group),paste0("MUTATIONS","_", figure, sep = "")))
  fileName = file_path_build(folder2Save,c(sampleDetail$Sample_ID,"MUTATIONS",figure),"bed")
  # mutations$SAMPLEID <- mutations$Sample_ID
  dumpSampleAsBedFile(
    dataToDump = mutations,
    fileName = fileName
  )

  # browser()
  folder2Save <- dir_check_and_create(resultFolderData,c(as.character(sampleDetail$Sample_Group),paste0("LESIONS","_", figure, sep = "")))
  fileName = file_path_build(folder2Save,c(sampleDetail$Sample_ID,"LESIONS",figure),"bed")
  lesions <- unique(lesions)
  if(nrow(lesions)>0)
  {
    # lesions$SAMPLEID <- sampleDetail$Sample_ID
    dumpSampleAsBedFile(
      dataToDump = lesions,
      fileName = fileName
    )
  }

  result["LESIONS_BOTH"] <- if (!is.null(lesions)) nrow(lesions) else 0
  result["MUTATIONS_BOTH"] <- if (!is.null(mutations)) nrow(mutations) else 0

  end_time_single_sample <- Sys.time()
  time_taken <- end_time_single_sample - start_time_single_sample
  message(sampleDetail$Sample_ID, " ", "Completed sample ", time_taken)

  return(result)
}
