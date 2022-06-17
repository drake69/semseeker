#' init environment
#'
#' @param result_folder where result of semseeker will bestored
#' @param parallel_strategy which strategy to use for parallel executio see future vignete: possibile values, none, multisession,sequential, multicore, cluster
#' @param maxResources percentage of how many available cores will be used default 90 percent, rounded to the lowest integer
#'
#' @return the working environment
init_env <- function(result_folder, maxResources = 90, parallel_strategy = "multisession", ...)
{
  set.seed(7658776)

  #allow export of object of 8gb with future
  options(future.globals.maxSize= 4 * 1024^3)

  ssEnv <- list()
  ssEnv$result_folderData <-  dir_check_and_create(result_folder, "Data")
  ssEnv$result_folderChart <-    dir_check_and_create(result_folder, "Chart")
  ssEnv$result_folderInference <-    dir_check_and_create(result_folder, "Inference")
  ssEnv$result_folderEuristic <-  dir_check_and_create(result_folder,"Euristic")
  random_file_name <- paste(stringi::stri_rand_strings(1, 7, pattern = "[A-Za-z0-9]"),".log", sep="")

  ssEnv$logFolder <-  dir_check_and_create(result_folder,c("log"))
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

  future::plan( future::sequential)
  if(parallel_strategy=="multisession")
    future::plan( future::multisession, workers = nCore)


  figures <- if(!exists("figures")) c("HYPO", "HYPER", "BOTH") else figures
  anomalies <- if(!exists("anomalies")) c("MUTATIONS","LESIONS","DELTAS") else anomalies
  metaareas <- if(!exists("metaareas")) c("GENE","ISLAND","DMR","CHR") else metaareas

  ssEnv$gene_subareas <- c("Body","TSS1500","5UTR","TSS200","1stExon","3UTR","ExonBnd","Whole")
  ssEnv$island_subareas <- c("N_Shore","S_Shore","N_Shelf","S_Shelf","Island", "Whole")

  ssEnv$keys_populations <-  data.frame("POPULATION"=c("Reference","Control","Case"))
  ssEnv$keys_figures <-  data.frame("FIGURE"=figures)
  ssEnv$keys_anomalies <-  data.frame("ANOMALY"=anomalies)
  ssEnv$keys_metaareas <- data.frame("METAAREA"=metaareas)

  ssEnv$keys_areas_island <-  expand.grid("GROUP"="ISLAND",
                                          "SUBGROUP"=island_subareas
  )
  ssEnv$keys_areas_gene <- expand.grid("GROUP"="GENE",
                                       "SUBGROUP"=gene_subareas
  )
  ssEnv$keys_areas_dmr <- expand.grid("GROUP"="DMR",
                                      "SUBGROUP"="DMR")

  ssEnv$keys_anomalies_figures_areas <- rbind(
    expand.grid("GROUP"="DMR",
                "SUBGROUP"="DMR",
                "FIGURE"=figures,
                "ANOMALY"=anomalies),
    expand.grid("GROUP"="GENE",
                "SUBGROUP"= gene_subareas,
                "FIGURE"=figures,
                "ANOMALY"=anomalies),
    expand.grid("GROUP"="ISLAND",
                "SUBGROUP"= island_subareas,
                "FIGURE"=figures,
                "ANOMALY"=anomalies)
  )

  ssEnv$keys <-  expand.grid("figures"=ssEnv$keys_figures[,1],"anomalies"=ssEnv$keys_anomalies[,1])

  probes_subGroups <- gene_subareas
  probes_Prefix <- "PROBES_Gene_"
  probes_MainGroupLabel <-  "GENE"
  probes_SubGroupLabel <- "GROUP"

  ssEnv$keys_gene_probes <-  expand.grid("prefix"=probes_Prefix,"maingrouplable"= probes_MainGroupLabel,"subgrouplable"= probes_SubGroupLabel,"subgroups"= probes_subGroups)
  # probes <-  expand.grid("prefix"=probes_Prefix,"maingrouplable"= probes_MainGroupLabel,"subgrouplable"= probes_SubGroupLabel,"subgroups"= probes_subGroups)

  # probes.27k
  # probes.450k
  # probes.850k

  probes_Prefix <- "PROBES_Island_"
  probes_subGroups <- island_subareas
  probes_MainGroupLabel <- "ISLAND"
  probes_SubGroupLabel <- "RELATION_TO_CPGISLAND"
  ssEnv$keys_island_probes <-  expand.grid("prefix"=probes_Prefix,"maingrouplable"= probes_MainGroupLabel,"subgrouplable"= probes_SubGroupLabel,"subgroups"= probes_subGroups)
  # probes <-  rbind(expand.grid("prefix"=probes_Prefix,"maingrouplable"= probes_MainGroupLabel,"subgrouplable"= probes_SubGroupLabel,"subgroups"= probes_subGroups), probes)

  probes_subGroups <- c("DMR")
  probes_Prefix <- "PROBES_DMR_"
  probes_MainGroupLabel <-  "DMR"
  probes_SubGroupLabel <- "GROUP"
  ssEnv$keys_dmr_probes <-  expand.grid("prefix"=probes_Prefix,"maingrouplable"= probes_MainGroupLabel,"subgrouplable"= probes_SubGroupLabel,"subgroups"= probes_subGroups)

  # probes <-  rbind( ssEnv$keys_island_probes, ssEnv$keys_gene_probes, ssEnv$keys_dmr_probes )

  ssEnv$functionToExport <- c( "analyze_single_sample",
                               "dump_sample_as_bed_file", "delta_single_sample","dir_check_and_create",
                               "file_path_build","analyze_single_sample_both",
                               "sort_by_chr_and_start", "test_match_order", "lesions_get",
                               "mutations_get","PROBES_Gene_3UTR", "PROBES_Gene_5UTR","PROBES_DMR_DMR","PROBES_Gene_Body")

  return(ssEnv)
}
