#' init environment
#'
#' @param resultFolder where result of semseeker will bestored
#' @param maxResources percentage of how many available cores will be used default 90 percent, rounded to the lowest integer
#' @return
#' @export
init_env <- function(resultFolder, maxResources = 90)
{

  if (.Platform$OS.type == "windows") {
    withAutoprint({
      utils::memory.size()
      utils::memory.size(TRUE)
      utils::memory.limit(16000)
    })
  }

  # setClass("employee", slots=list(name="character", id="numeric", contact="character"))
  # obj <- new("employee",name="Steven", id=1002, contact="West Avenue")

  resultFolderData <<- dir_check_and_create(resultFolder, "Data")
  resultFolderChart <<-   dir_check_and_create(resultFolder, "Chart")
  resultFolderInference <<-   dir_check_and_create(resultFolder, "Inference")
  resultFolderEuristic <<- dir_check_and_create(resultFolder,"Euristic")
  folderLog <<- dir_check_and_create("/tmp",c("semseeker","log"))

  # bootstrap cluster
  nCore <- parallel::detectCores(all.tests = FALSE, logical = TRUE)
  nCore <- floor(nCore * maxResources/100 )
  outFile <- file.path(folderLog, "cluster_r.out")
  # print(outFile)
  computationCluster <<- parallel::makeCluster(parallel::detectCores(all.tests = FALSE, logical = TRUE) - 1, type = "PSOCK", outfile = outFile)
  doParallel::registerDoParallel(computationCluster)


  keys.populations <<- c("Reference","Control","Case")
  keys.figures <<- c("HYPO", "HYPER", "BOTH")
  keys.anomalies <<- c("MUTATIONS","LESIONS")

  keys <<- expand.grid("figures"=keys.figures,"anomalies"=keys.anomalies)

  probes.subGroups <- c("Body","TSS1500","5UTR","TSS200","1stExon","3UTR","ExonBnd","Whole")
  probes.Prefix = "PROBES_Gene_"
  probes.MainGroupLabel =  "GENE"
  probes.SubGroupLabel="GROUP"
  probes <<- expand.grid("prefix"=probes.Prefix,"maingrouplable"= probes.MainGroupLabel,"subgrouplable"= probes.SubGroupLabel,"subgroups"= probes.subGroups)

  # probes.27k
  # probes.450k
  # probes.850k

  probes.Prefix <- "PROBES_Island_"
  probes.subGroups <- c("N_Shore","S_Shore","N_Shelf","S_Shelf","Island", "Whole")
  probes.MainGroupLabel <- "ISLAND"
  probes.SubGroupLabel <- "RELATION_TO_CPGISLAND"
  probes <<- rbind(expand.grid("prefix"=probes.Prefix,"maingrouplable"= probes.MainGroupLabel,"subgrouplable"= probes.SubGroupLabel,"subgroups"= probes.subGroups), probes)

  probes.subGroups <- c("DMR")
  probes.Prefix = "PROBES_DMR_"
  probes.MainGroupLabel =  "DMR"
  probes.SubGroupLabel="GROUP"
  probes <<- rbind(expand.grid("prefix"=probes.Prefix,"maingrouplable"= probes.MainGroupLabel,"subgrouplable"= probes.SubGroupLabel,"subgroups"= probes.subGroups), probes)

}
