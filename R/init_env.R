#' init environment
#' @param result_folder where result of semseeker will bestored
#' @param maxResources percentage of how many available cores will be used default 90 percent, rounded to the lowest integer
#' @return the working environment
init_env <- function(result_folder, maxResources = 90)
{


  # if(getRversion() >= "2.15.1")  utils::globalVariables(c( "chartFolder","ssEnv$logFolder","computationCluster",
  #                                                          "ssEnv$result_folderChart","ssEnv$result_folderInference","ssEnv$result_folderEuristic",
  #                                                          "ssEnv$keys_anomalies","ssEnv$keys_figures","ssEnv$keys_populations","ssEnv$keys",
  #                                                          "probes","functionToExport"), add = FALSE)

  # if(length(utils::globalVariables())!=0)
  #   return()

  # if(getRversion() >= "2.15.1")
  #   utils::globalVariables( names= c( "ssEnv$result_folderData","chartFolder","ssEnv$logFolder","computationCluster",
  #                             "ssEnv$result_folderChart","ssEnv$result_folderInference","ssEnv$result_folderEuristic",
  #                             "ssEnv$keys_anomalies","ssEnv$keys_figures","ssEnv$keys_populations","ssEnv$keys",
  #                             "probes"), add = FALSE, package = "semseeker")


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

  # assign("ssEnv$result_folderData", dir_check_and_create(result_folder, "Data"))

  ssEnv <- list()
  ssEnv$result_folderData <-  dir_check_and_create(result_folder, "Data")
  ssEnv$result_folderChart <-    dir_check_and_create(result_folder, "Chart")
  ssEnv$result_folderInference <-    dir_check_and_create(result_folder, "Inference")
  ssEnv$result_folderEuristic <-  dir_check_and_create(result_folder,"Euristic")
  ssEnv$logFolder <-  dir_check_and_create("/tmp",c("semseeker","log"))
  foreachIndex <- 0

  # print(outFile)
  # computationCluster <- parallel::makeCluster(nCore, type = "PSOCK", outfile = outFile)

  # doFuture::registerDoFuture()
  # computationCluster <- doFuture::makeCluster(nCore, outFile = outfile)
  # plan(cluster, workers = computationCluster)
  chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")

  if (nzchar(chk) && chk == "TRUE") {
    # use 2 cores in CRAN/Travis/AppVeyor
    nCore <- 2L
  } else {
    # use all cores in devtools::test()
    nCore <- future::availableCores()
    nCore <- floor(nCore * maxResources/100 )
  }
  # bootstrap cluster
  outFile <- file.path(ssEnv$logFolder, "cluster_r.out")

  options(doFuture.foreach.export = ".export-and-automatic-with-warning")
  doFuture::registerDoFuture()
  future::plan( future::multisession, workers = nCore)


  ssEnv$keys_populations <-  c("Reference","Control","Case")
  ssEnv$keys_figures <-  c("HYPO", "HYPER", "BOTH")
  ssEnv$keys_anomalies <-  c("MUTATIONS","LESIONS")

  ssEnv$keys <-  expand.grid("figures"=ssEnv$keys_figures,"anomalies"=ssEnv$keys_anomalies)

  probes_subGroups <- c("Body","TSS1500","5UTR","TSS200","1stExon","3UTR","ExonBnd","Whole")
  probes_Prefix <- "PROBES_Gene_"
  probes_MainGroupLabel <-  "GENE"
  probes_SubGroupLabel <- "GROUP"

  ssEnv$keys_gene_probes <-  expand.grid("prefix"=probes_Prefix,"maingrouplable"= probes_MainGroupLabel,"subgrouplable"= probes_SubGroupLabel,"subgroups"= probes_subGroups)
  # probes <-  expand.grid("prefix"=probes_Prefix,"maingrouplable"= probes_MainGroupLabel,"subgrouplable"= probes_SubGroupLabel,"subgroups"= probes_subGroups)

  # probes.27k
  # probes.450k
  # probes.850k

  probes_Prefix <- "PROBES_Island_"
  probes_subGroups <- c("N_Shore","S_Shore","N_Shelf","S_Shelf","Island", "Whole")
  probes_MainGroupLabel <- "ISLAND"
  probes_SubGroupLabel <- "RELATION_TO_CPGISLAND"
  ssEnv$keys_island_probes <-  expand.grid("prefix"=probes_Prefix,"maingrouplable"= probes_MainGroupLabel,"subgrouplable"= probes_SubGroupLabel,"subgroups"= probes_subGroups)
  # probes <-  rbind(expand.grid("prefix"=probes_Prefix,"maingrouplable"= probes_MainGroupLabel,"subgrouplable"= probes_SubGroupLabel,"subgroups"= probes_subGroups), probes)

  probes_subGroups <- c("DMR")
  probes_Prefix <- "PROBES_DMR_"
  probes_MainGroupLabel <-  "DMR"
  probes_SubGroupLabel <- "GROUP"
  ssEnv$keys_dmr_probes <-  expand.grid("prefix"=probes_Prefix,"maingrouplable"= probes_MainGroupLabel,"subgrouplable"= probes_SubGroupLabel,"subgroups"= probes_subGroups)

  probes <-  rbind( ssEnv$keys_island_probes, ssEnv$keys_gene_probes, ssEnv$keys_dmr_probes )


  # parallel::clusterExport(envir=my_env, cl=computationCluster, varlist = list("ssEnv$result_folderData","ssEnv$logFolder","computationCluster",
  #                                         "ssEnv$result_folderChart","ssEnv$result_folderInference","ssEnv$result_folderEuristic",
  #                                         "ssEnv$keys_anomalies","ssEnv$keys_figures","ssEnv$keys_populations","ssEnv$keys",
  #                                         "probes"))
  #
  # parallel::clusterExport(envir=environment(),
  #                         cl=computationCluster,
  #                         varlist = list( "analyze_single_sample",
  #                           "dump_sample_as_bed_file", "delta_single_sample","dir_check_and_create",
  #                           "file_path_build","analyze_single_sample_both",
  #                           "createPivotResultFromMultipleBed", "sort_by_chr_and_start", "test_match_order", "lesions_get",
  #                           "mutations_get","PROBES_Gene_3UTR", "PROBES_Gene_5UTR","PROBES_DMR_DMR","PROBES_Gene_Body",
  #                           "ssEnv$result_folderData","ssEnv$logFolder","computationCluster",
  #                           "ssEnv$result_folderChart","ssEnv$result_folderInference","ssEnv$result_folderEuristic",
  #                           "ssEnv$keys_anomalies","ssEnv$keys_figures","ssEnv$keys_populations","ssEnv$keys",
  #                           "probes"))

  # ,
  # "ssEnv$result_folderData","ssEnv$logFolder","computationCluster",
  # "ssEnv$result_folderChart","ssEnv$result_folderInference","ssEnv$result_folderEuristic",
  # "ssEnv$keys_anomalies","ssEnv$keys_figures","ssEnv$keys_populations","ssEnv$keys",
  # "probes"

  ssEnv$functionToExport <- c( "analyze_single_sample",
                            "dump_sample_as_bed_file", "delta_single_sample","dir_check_and_create",
                            "file_path_build","analyze_single_sample_both",
                            "sort_by_chr_and_start", "test_match_order", "lesions_get",
                            "mutations_get","PROBES_Gene_3UTR", "PROBES_Gene_5UTR","PROBES_DMR_DMR","PROBES_Gene_Body")

  # options(doFuture.foreach.export = ".export-and-automatic-with-warning")

  #
  #
  # parallel::clusterExport(envir=my_env, cl=computationCluster, varlist = ls(my_env))

  # lockBinding("ssEnv$result_folderData","chartFolder","ssEnv$logFolder","computationCluster",
  #             "ssEnv$result_folderChart","ssEnv$result_folderInference","ssEnv$result_folderEuristic",
  #             "ssEnv$keys_anomalies","ssEnv$keys_figures","ssEnv$keys_populations","ssEnv$keys",
  #             "probes", e)

  return(ssEnv)
}
