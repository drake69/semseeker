#' init environment
#' @param resultFolder where result of semseeker will bestored
#' @param maxResources percentage of how many available cores will be used default 90 percent, rounded to the lowest integer
#' @return null
init_env <- function(resultFolder, maxResources = 90)
{

  # if(length(utils::globalVariables())!=0)
  #   return()

  if(getRversion() >= "2.15.1")
    utils::globalVariables( names= c( "resultFolderData","chartFolder","logFolder","computationCluster",
                              "resultFolderChart","resultFolderInference","resultFolderEuristic",
                              "keys_anomalies","keys_figures","keys_populations","keys",
                              "probes"), add = FALSE, package = "semseeker")


  if (.Platform$OS.type == "windows") {
    withAutoprint({
      utils::memory.size()
      utils::memory.size(TRUE)
      utils::memory.limit(16000)
    })
  }

  # setClass("employee", slots=list(name="character", id="numeric", contact="character"))
  # obj <- new("employee",name="Steven", id=1002, contact="West Avenue")

  # e <- new.env()

  resultFolderData  <<-  dir_check_and_create(resultFolder, "Data")
  resultFolderChart  <<-    dir_check_and_create(resultFolder, "Chart")
  resultFolderInference  <<-    dir_check_and_create(resultFolder, "Inference")
  resultFolderEuristic  <<-  dir_check_and_create(resultFolder,"Euristic")
  logFolder  <<-  dir_check_and_create("/tmp",c("semseeker","log"))

  # bootstrap cluster
  nCore <- parallel::detectCores(all.tests = FALSE, logical = TRUE)
  nCore <- floor(nCore * maxResources/100 )
  outFile <- file.path(logFolder, "cluster_r.out")
  # print(outFile)
  computationCluster  <<-  parallel::makeCluster(parallel::detectCores(all.tests = FALSE, logical = TRUE) - 1, type = "PSOCK", outfile = outFile)
  doParallel::registerDoParallel(computationCluster)


  keys_populations  <<-  c("Reference","Control","Case")
  keys_figures  <<-  c("HYPO", "HYPER", "BOTH")
  keys_anomalies  <<-  c("MUTATIONS","LESIONS")

  keys  <<-  expand.grid("figures"=keys_figures,"anomalies"=keys_anomalies)

  probes.subGroups <- c("Body","TSS1500","5UTR","TSS200","1stExon","3UTR","ExonBnd","Whole")
  probes.Prefix <- "PROBES_Gene_"
  probes.MainGroupLabel <-  "GENE"
  probes.SubGroupLabel <- "GROUP"

  keys_gene_probes  <-  expand.grid("prefix"=probes.Prefix,"maingrouplable"= probes.MainGroupLabel,"subgrouplable"= probes.SubGroupLabel,"subgroups"= probes.subGroups)
  # probes  <<-  expand.grid("prefix"=probes.Prefix,"maingrouplable"= probes.MainGroupLabel,"subgrouplable"= probes.SubGroupLabel,"subgroups"= probes.subGroups)

  # probes.27k
  # probes.450k
  # probes.850k

  probes.Prefix <- "PROBES_Island_"
  probes.subGroups <- c("N_Shore","S_Shore","N_Shelf","S_Shelf","Island", "Whole")
  probes.MainGroupLabel <- "ISLAND"
  probes.SubGroupLabel <- "RELATION_TO_CPGISLAND"
  keys_island_probes  <-  expand.grid("prefix"=probes.Prefix,"maingrouplable"= probes.MainGroupLabel,"subgrouplable"= probes.SubGroupLabel,"subgroups"= probes.subGroups)
  # probes  <<-  rbind(expand.grid("prefix"=probes.Prefix,"maingrouplable"= probes.MainGroupLabel,"subgrouplable"= probes.SubGroupLabel,"subgroups"= probes.subGroups), probes)

  probes.subGroups <- c("DMR")
  probes.Prefix <- "PROBES_DMR_"
  probes.MainGroupLabel <-  "DMR"
  probes.SubGroupLabel <- "GROUP"
  keys_dmr_probes  <-  expand.grid("prefix"=probes.Prefix,"maingrouplable"= probes.MainGroupLabel,"subgrouplable"= probes.SubGroupLabel,"subgroups"= probes.subGroups)
  # probes  <<-  rbind(expand.grid("prefix"=probes.Prefix,"maingrouplable"= probes.MainGroupLabel,"subgrouplable"= probes.SubGroupLabel,"subgroups"= probes.subGroups), probes)


  parallel::clusterExport(envir=environment(), cl=computationCluster, varlist = list( "analyzeSingleSample", "dumpSampleAsBedFile", "deltaSingleSample","dir_check_and_create","resultFolderData",
                                          "file_path_build","analyzeSingleSampleBoth",
                                          "createPivotResultFromMultipleBed", "sortByCHRandSTART", "test_match_order", "getLesions",
                                          "getMutations","PROBES_Gene_3UTR", "PROBES_Gene_5UTR","PROBES_DMR_DMR","PROBES_Gene_Body",
                                          "resultFolderData","logFolder","computationCluster",
                                          "resultFolderChart","resultFolderInference","resultFolderEuristic",
                                          "keys_anomalies","keys_figures","keys_populations","keys",
                                          "probes"))

  # lockBinding("resultFolderData","chartFolder","logFolder","computationCluster",
  #             "resultFolderChart","resultFolderInference","resultFolderEuristic",
  #             "keys_anomalies","keys_figures","keys_populations","keys",
  #             "probes", e)

}
