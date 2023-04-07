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

  # utils::data("PROBES")
  # utils::data("PROBES_CHR_CHR")
  set.seed(7658776)

  #allow export of object of 8gb with future
  options(future.globals.maxSize= 16 * 1024^3)

  ssEnv <- list()
  ssEnv$parallel_strategy <- parallel_strategy
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
    nCore <- if(floor(future::availableCores() * maxResources/100 ) > nCore ) nCore else floor(future::availableCores() * maxResources/100 )
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
      message ("ERROR: ", Sys.time(), " Function's arguments must be passed explicitily !")
      message(cond)
      stop()
    }
  )
  arguments <- list(...)

  ssEnv$showprogress <- FALSE
  if(!is.null(arguments[["showprogress"]]))
  {
    ssEnv$showprogress <- if(is.null(arguments[["showprogress"]])) ssEnv$keys_figures_default[,1] else arguments$showprogress
  }

  # TODO: improve planning parallel management using also cluster
  if(parallel_strategy=="multisession")
  {
    future::plan( future::multisession, workers = nCore)
    message("INFO: ", Sys.time(), " I will work in multisession with:", nCore, " Cores")
  }
  if(parallel_strategy=="multicore")
  {
    future::plan( future::multicore, workers = nCore)
    message("INFO: ", Sys.time(), " I will work in muticore with:", nCore," Cores")
  }
  if(parallel_strategy=="cluster")
  {
    message ("ERROR: ", Sys.time(), " Cluster feature not implemented!")
    stop()
    future::plan( future::cluster, workers = nCore)
  }
  if(parallel_strategy!="multisession" & parallel_strategy!="multicore"
     & parallel_strategy!="cluster")
  {
    future::plan( future::sequential)
    message("INFO: ", Sys.time(), " I will work in sequential mode")
  }

  ssEnv$keys_figures_default <-  data.frame("FIGURE"=c("HYPO", "HYPER", "BOTH"))
  ssEnv$keys_anomalies_default <-  data.frame("ANOMALY"=c("MUTATIONS","LESIONS","DELTAS","DELTAQ"))
  ssEnv$keys_metaareas_default <- data.frame("METAAREA"=c("GENE","ISLAND","DMR","CHR","PROBE"))
  # ,"PROBE"

  figures <- if(is.null(arguments[["figures"]])) ssEnv$keys_figures_default[,1] else arguments$figures
  anomalies <- if(is.null(arguments[["anomalies"]]))  ssEnv$keys_anomalies_default[,1] else arguments$anomalies
  metaareas <- if(is.null(arguments[["metaareas"]]))  ssEnv$keys_metaareas_default[,1] else arguments$metaareas

  if(sum(figures %in% ssEnv$keys_figures_default[,1])==0)
  {
    message("INFO: ", Sys.time(), " The only allowed figures values are:", ssEnv$keys_figures_default)
    stop()
  }
  if(sum(anomalies %in% ssEnv$keys_anomalies_default[,1])==0)
  {
    message("INFO: ", Sys.time(), " The only allowed anomalies values are:", ssEnv$keys_anomalies_default)
    stop()
  }
  if(sum(metaareas %in% ssEnv$keys_metaareas_default[,1])==0)
  {
    message("INFO: ", Sys.time(), " The only allowed areas values are:", ssEnv$keys_metaareas_default)
    stop()
  }

  message("INFO: ", Sys.time(), " I will focus on:", paste(anomalies, collapse = " ", sep =" "), " due to ",  paste(figures, collapse = " ", sep =" "), " of ",  paste(metaareas, collapse = " ", sep =" "))

  ssEnv$gene_subareas <- data.frame("subarea"=c("Body","TSS1500","5UTR","TSS200","1stExon","3UTR","5UTR","ExonBnd","Whole"))
  ssEnv$island_subareas <- data.frame("subarea"=c("N_Shore","S_Shore","N_Shelf","S_Shelf","Island", "Whole"))

  if(!is.null(arguments[["metaareassub"]]))
  {
    ssEnv$gene_subareas <- as.data.frame(ssEnv$gene_subareas[ssEnv$gene_subareas[,"subarea"]  %in% arguments$metaareassub,])
    ssEnv$island_subareas <- as.data.frame(ssEnv$island_subareas[ssEnv$island_subareas[,"subarea"] %in% arguments$metaareassub,  ])
  }

  ssEnv$keys_populations <-  data.frame("POPULATION"=c("Reference","Control","Case"))
  ssEnv$keys_figures <-  data.frame("FIGURE"=figures)
  ssEnv$keys_anomalies <-  data.frame("ANOMALY"=anomalies)
  ssEnv$keys_metaareas <- data.frame("METAAREA"=metaareas)

  if("ISLAND" %in% metaareas)
    ssEnv$keys_areas_island <-  expand.grid("GROUP"="ISLAND",
                                            "SUBGROUP"=ssEnv$island_subareas[,1]
    )
  else
  {
    ssEnv$keys_areas_island <-  expand.grid("GROUP"="ISLAND","SUBGROUP"="")[-1,]
  }

  if("GENE" %in% metaareas)
    ssEnv$keys_areas_gene <- expand.grid("GROUP"="GENE",
                                       "SUBGROUP"=ssEnv$gene_subareas[,1]
    )
  else
  {
    ssEnv$keys_areas_gene <- expand.grid("GROUP"="GENE","SUBGROUP"="")[-1,]
  }

  ssEnv$keys_areas_dmr <- expand.grid("GROUP"="DMR","SUBGROUP"="DMR")
  if (!("DMR" %in% metaareas))
    ssEnv$keys_areas_dmr <- ssEnv$keys_areas_dmr[-1,]

  ssEnv$keys_areas_chr <- expand.grid("GROUP"="CHR","SUBGROUP"="CHR")
  if (!("CHR" %in% metaareas))
    ssEnv$keys_areas_chr <- ssEnv$keys_areas_chr[-1,]

  ssEnv$keys_areas_probe <- expand.grid("GROUP"="PROBE","SUBGROUP"="PROBE")
  if(!("PROBE" %in% metaareas))
    ssEnv$keys_areas_probe <-  ssEnv$keys_areas_probe[-1,]

  ssEnv$keys_anomalies_figures_areas <- rbind(
    expand.grid("GROUP"="CHR",
                "SUBGROUP"="CHR",
                "FIGURE"=figures,
                "ANOMALY"=anomalies),
    expand.grid("GROUP"="PROBE",
                "SUBGROUP"="PROBE",
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

  # remove anomaies if metaareas not in options
  ssEnv$keys_anomalies_figures_areas <- ssEnv$keys_anomalies_figures_areas[ ssEnv$keys_anomalies_figures_areas$GROUP %in% metaareas, ]

  ssEnv$keys <-  expand.grid("figures"=ssEnv$keys_figures[,1],"anomalies"=ssEnv$keys_anomalies[,1])

  if(nrow(ssEnv$keys_areas_probe)>0)
  {
    probes_subGroups <- ""
    probes_Prefix <- "PROBES"
    probes_MainGroupLabel <-  "PROBE"
    probes_SubGroupLabel <- "GROUP"
  }
  else
  {
    probes_subGroups <- NULL
    probes_Prefix <- NULL
    probes_MainGroupLabel <-  NULL
    probes_SubGroupLabel <- NULL
  }

  ssEnv$keys_probe_probes <-  expand.grid("prefix"=probes_Prefix,"maingrouplable"= probes_MainGroupLabel,"subgrouplable"= probes_SubGroupLabel,"subgroups"= probes_subGroups)

  if(nrow(ssEnv$keys_areas_chr)>0)
  {
    probes_subGroups <- "CHR"
    probes_Prefix <- "PROBES_CHR_"
    probes_MainGroupLabel <-  "CHR"
    probes_SubGroupLabel <- "GROUP"
  }
  else
  {
    probes_subGroups <- NULL
    probes_Prefix <- NULL
    probes_MainGroupLabel <-  NULL
    probes_SubGroupLabel <- NULL
  }

  ssEnv$keys_chr_probes <-  expand.grid("prefix"=probes_Prefix,"maingrouplable"= probes_MainGroupLabel,"subgrouplable"= probes_SubGroupLabel,"subgroups"= probes_subGroups)

  if(nrow(ssEnv$keys_areas_gene)>0)
  {
    probes_subGroups <- ssEnv$gene_subareas[,1]
    probes_Prefix <- "PROBES_Gene_"
    probes_MainGroupLabel <-  "GENE"
    probes_SubGroupLabel <- "GROUP"
  }
  else
  {
    probes_subGroups <- NULL
    probes_Prefix <- NULL
    probes_MainGroupLabel <-  NULL
    probes_SubGroupLabel <- NULL
  }

  ssEnv$keys_gene_probes <-  expand.grid("prefix"=probes_Prefix,"maingrouplable"= probes_MainGroupLabel,"subgrouplable"= probes_SubGroupLabel,"subgroups"= probes_subGroups)
  # probes <-  expand.grid("prefix"=probes_Prefix,"maingrouplable"= probes_MainGroupLabel,"subgrouplable"= probes_SubGroupLabel,"subgroups"= probes_subGroups)

  # probes.27k
  # probes.450k
  # probes.850k

  if(nrow(ssEnv$keys_areas_island)>0)
  {
    probes_Prefix <- "PROBES_Island_"
    probes_subGroups <- ssEnv$island_subareas[,1]
    probes_MainGroupLabel <- "ISLAND"
    probes_SubGroupLabel <- "RELATION_TO_CPGISLAND"
  }
  else
  {
    probes_subGroups <- NULL
    probes_Prefix <- NULL
    probes_MainGroupLabel <-  NULL
    probes_SubGroupLabel <- NULL
  }
  ssEnv$keys_island_probes <-  expand.grid("prefix"=probes_Prefix,"maingrouplable"= probes_MainGroupLabel,"subgrouplable"= probes_SubGroupLabel,"subgroups"= probes_subGroups)
  # probes <-  rbind(expand.grid("prefix"=probes_Prefix,"maingrouplable"= probes_MainGroupLabel,"subgrouplable"= probes_SubGroupLabel,"subgroups"= probes_subGroups), probes)

  if(nrow(ssEnv$keys_areas_dmr)>0)
  {
    probes_subGroups <- c("DMR")
    probes_Prefix <- "PROBES_DMR_"
    probes_MainGroupLabel <-  "DMR"
    probes_SubGroupLabel <- "GROUP"
  }
  else
  {
    probes_subGroups <- NULL
    probes_Prefix <- NULL
    probes_MainGroupLabel <-  NULL
    probes_SubGroupLabel <- NULL
  }
  ssEnv$keys_dmr_probes <-  expand.grid("prefix"=probes_Prefix,"maingrouplable"= probes_MainGroupLabel,"subgrouplable"= probes_SubGroupLabel,"subgroups"= probes_subGroups)

  # probes <-  rbind( ssEnv$keys_island_probes, ssEnv$keys_gene_probes, ssEnv$keys_dmr_probes )

  ssEnv$functionToExport <- c( "analyze_single_sample",
                               "dump_sample_as_bed_file", "delta_single_sample","dir_check_and_create",
                               "file_path_build","analyze_single_sample_both",
                               "sort_by_chr_and_start", "test_match_order", "lesions_get",
                               "mutations_get","PROBES_Gene_3UTR", "PROBES_Gene_5UTR","PROBES_DMR_DMR","PROBES_Gene_Body")


  # to manage progress bar
  if(ssEnv$showprogress)
  {
    progressr::handlers(global = TRUE)
    progressr::handlers("progress")
  }
  return(ssEnv)
}
