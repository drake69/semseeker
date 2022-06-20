#' init environment
#'
#' @param result_folder where result of semseeker will bestored
#' @param parallel_strategy which strategy to use for parallel executio see future vignete: possibile values, none, multisession,sequential, multicore, cluster
#' @param maxResources percentage of how many available cores will be used default 90 percent, rounded to the lowest integer
#' @param ... other options to filter elaborations
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

  chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")

  if (nzchar(chk) && chk == "TRUE") {
    # use 2 cores in CRAN/Travis/AppVeyor
    nCore <- 2L
  } else {
    # use all cores in devtools::test()
    nCore <- future::availableCores() - 1
    nCore <- floor(nCore * maxResources/100 )
  }
  # bootstrap cluster
  outFile <- file.path(ssEnv$logFolder, "cluster_r.out")

  options(doFuture.foreach.export = ".export-and-automatic-with-warning")
  doFuture::registerDoFuture()


  tryCatch(
    {
      test_it <- list(...)
    },
    error = function(cond)  {
      message ("Function's arguments must be passed explicitily !")
      message(cond)
      stop()
    }
  )
  arguments <- list(...)

  # TODO: improve planning parallel management using also cluster
  future::plan( future::sequential)
  if(parallel_strategy=="multisession")
  {
    future::plan( future::multisession, workers = nCore)
    message("I will work in multisession with:", nCore, " Cores")
  }
  if(parallel_strategy=="multicore")
  {
    future::plan( future::multicore, workers = nCore," Cores")
    message("I will work in muticore with:", nCore)
  }
  if(parallel_strategy=="cluster")
  {
    message ("Cluster feature not implemented!")
    stop()
    future::plan( future::cluster, workers = nCore)
  }
  if(parallel_strategy!="multisession" & parallel_strategy!="multicore"
     & parallel_strategy!="cluster")
  {
    future::plan( future::multisession, workers = nCore)
    message("I will work in sequential mode")
  }

  ssEnv$keys_figures_default <-  data.frame("FIGURE"=c("HYPO", "HYPER", "BOTH"))
  ssEnv$keys_anomalies_default <-  data.frame("ANOMALY"=c("MUTATIONS","LESIONS","DELTAS","DELTAQ"))
  ssEnv$keys_metaareas_default <- data.frame("METAAREA"=c("GENE","ISLAND","DMR","CHR"))



  figures <- if(is.null(arguments[["figures"]])) ssEnv$keys_figures_default[,1] else arguments$figures
  anomalies <- if(is.null(arguments[["anomalies"]]))  ssEnv$keys_anomalies_default[,1] else arguments$anomalies
  metaareas <- if(is.null(arguments[["metaareas"]]))  ssEnv$keys_metaareas_default[,1] else arguments$metaareas

  if(sum(figures %in% ssEnv$keys_figures_default[,1])==0)
  {
    message("The only allowed figures values are:", ssEnv$keys_figures_default)
    stop()
  }
  if(sum(anomalies %in% ssEnv$keys_anomalies_default[,1])==0)
  {
    message("The only allowed anomalies values are:", ssEnv$keys_anomalies_default)
    stop()
  }
  if(sum(metaareas %in% ssEnv$keys_metaareas_default[,1])==0)
  {
    message("The only allowed areas values are:", ssEnv$keys_metaareas_default)
    stop()
  }

  message("I will focus on:", paste(anomalies, collapse = " ", sep =" "), " due to ",  paste(figures, collapse = " ", sep =" "), " of ",  paste(metaareas, collapse = " ", sep =" "))

  ssEnv$gene_subareas <- data.frame(c("Body","TSS1500","5UTR","TSS200","1stExon","3UTR","ExonBnd","Whole"))
  ssEnv$island_subareas <- data.frame(c("N_Shore","S_Shore","N_Shelf","S_Shelf","Island", "Whole"))

  ssEnv$keys_populations <-  data.frame("POPULATION"=c("Reference","Control","Case"))
  ssEnv$keys_figures <-  data.frame("FIGURE"=figures)
  ssEnv$keys_anomalies <-  data.frame("ANOMALY"=anomalies)
  ssEnv$keys_metaareas <- data.frame("METAAREA"=metaareas)

  ssEnv$keys_areas_island <-  expand.grid("GROUP"="ISLAND",
                                          "SUBGROUP"=ssEnv$island_subareas[,1]
  )
  ssEnv$keys_areas_gene <- expand.grid("GROUP"="GENE",
                                       "SUBGROUP"=ssEnv$gene_subareas[,1]
  )
  ssEnv$keys_areas_dmr <- expand.grid("GROUP"="DMR",
                                      "SUBGROUP"="DMR")

  ssEnv$keys_areas_chr <- expand.grid("GROUP"="CHR",
                                      "SUBGROUP"="CHR")

  ssEnv$keys_anomalies_figures_areas <- rbind(
    expand.grid("GROUP"="CHR",
                "SUBGROUP"="CHR",
                "FIGURE"=figures,
                "ANOMALY"=anomalies),
    expand.grid("GROUP"="DMR",
                "SUBGROUP"="DMR",
                "FIGURE"=figures,
                "ANOMALY"=anomalies),
    expand.grid("GROUP"="GENE",
                "SUBGROUP"= ssEnv$gene_subareas[,1],
                "FIGURE"=figures,
                "ANOMALY"=anomalies),
    expand.grid("GROUP"="ISLAND",
                "SUBGROUP"= ssEnv$island_subareas[,1],
                "FIGURE"=figures,
                "ANOMALY"=anomalies)
  )

  ssEnv$keys_anomalies_figures_areas <- ssEnv$keys_anomalies_figures_areas[ ssEnv$keys_anomalies_figures_areas$GROUP %in% metaareas, ]

  ssEnv$keys <-  expand.grid("figures"=ssEnv$keys_figures[,1],"anomalies"=ssEnv$keys_anomalies[,1])

  probes_subGroups <- ssEnv$gene_subareas[,1]
  probes_Prefix <- "PROBES_Gene_"
  probes_MainGroupLabel <-  "GENE"
  probes_SubGroupLabel <- "GROUP"

  ssEnv$keys_gene_probes <-  expand.grid("prefix"=probes_Prefix,"maingrouplable"= probes_MainGroupLabel,"subgrouplable"= probes_SubGroupLabel,"subgroups"= probes_subGroups)
  # probes <-  expand.grid("prefix"=probes_Prefix,"maingrouplable"= probes_MainGroupLabel,"subgrouplable"= probes_SubGroupLabel,"subgroups"= probes_subGroups)

  # probes.27k
  # probes.450k
  # probes.850k

  probes_Prefix <- "PROBES_Island_"
  probes_subGroups <- ssEnv$island_subareas[,1]
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
