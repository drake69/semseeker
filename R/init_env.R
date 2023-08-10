#' init ssEnvonment
#'
#' @param result_folder where result of semseeker will bestored
#' @param parallel_strategy which strategy to use for parallel executio see future vignete: possibile values, none, multisession,sequential, multicore, cluster
#' @param maxResources percentage of how many available cores will be used default 90 percent, rounded to the lowest integer
#' @param ... other options to filter elaborations
#'
#' @return the working ssEnvonment
init_env <- function(result_folder, maxResources = 90, parallel_strategy = "multicore", ...)
{

  # utils::data("PROBES")
  # utils::data("PROBES_CHR_CHR")
  set.seed(7658776)

  #allow export of object of 32gb with future
  options(future.globals.maxSize= 32 * 1024^3)

  ssEnv <- get_session_info(result_folder)
  ssEnv$parallel_strategy <- parallel_strategy

  tmp <- tempdir()
  ssEnv$temp_folder <-  paste(tmp,"/semseeker/",stringi::stri_rand_strings(1, 7, pattern = "[A-Za-z0-9]"),sep="")
  ssEnv$result_folderData <-  dir_check_and_create(result_folder, "Data")
  ssEnv$result_folderChart <-    dir_check_and_create(result_folder, "Chart")
  ssEnv$result_folderInference <-    dir_check_and_create(result_folder, "Inference")
  ssEnv$result_folderEuristic <-  dir_check_and_create(result_folder,"Euristic")
  ssEnv$session_folder <-  dir_check_and_create(result_folder,c("Log"))
  random_file_name <- paste(stringi::stri_rand_strings(1, 7, pattern = "[A-Za-z0-9]"),".log", sep="")

  if (sink.number() != 0)
    sink(NULL)
  sink(file.path(ssEnv$session_folder,"session_output.log"), split = TRUE)

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
  outFile <- file.path(ssEnv$session_folder, "cluster_r.out")

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
    ssEnv$showprogress <- if(is.null(arguments[["showprogress"]])) keys_figures_default[,1] else arguments$showprogress
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
    stop("I'm STOPPING HERE!")
    future::plan( future::cluster, workers = nCore)
  }
  if(parallel_strategy!="multisession" & parallel_strategy!="multicore"
     & parallel_strategy!="cluster")
  {
    future::plan( future::sequential)
    message("INFO: ", Sys.time(), " I will work in sequential mode")
  }

  ssEnvTemp <- get_session_info(result_folder)
  if (!is.null(ssEnvTemp) & !length(ssEnvTemp)<2)
  {
    # we had to return aold session info this init could be calles by other than semseeker function, it means
    # we don't know which marker or figure were calculated
    message("Reusing old session info !")
    message("INFO: ", Sys.time(), " I will focus on: ", paste(unique(ssEnvTemp$keys_markers_figures[,"MARKER"]), collapse = " ", sep =" "), " due to ",
      paste(unique(ssEnvTemp$keys_markers_figures[,"FIGURE"]), collapse = " ", sep =" "),
      " of ",  paste(unique(ssEnvTemp$keys_areas_subareas[,"AREA"]), collapse = " ", sep =" "),
      " for ",  paste(unique(ssEnvTemp$keys_areas_subareas[,"SUBAREA"]), collapse = " ", sep =" "))
    if(dir.exists(ssEnv$session_folder))
    {
      update_session_info(ssEnvTemp)
      return(ssEnvTemp)
    }
  }

  # set default values
  keys_figures_default <-  data.frame("FIGURE"=c("HYPO", "HYPER", "BOTH"))
  keys_markers_default <-  data.frame("MARKER"=c("MUTATIONS","LESIONS","DELTAS","DELTAQ","DELTAR","DELTARQ","BETA"))
  keys_markers_default$EXT <-  c("bed","bed","bedgraph","bed","bedgraph","bed","bedgraph")
  keys_areas_default <- data.frame("areas"=c("GENE","ISLAND","DMR","CHR","PROBE"))
  keys_gene_subareas_default <- data.frame("subarea"=c("BODY","TSS1500","TSS200","1STEXON","3UTR","5UTR","EXONBND","WHOLE"))
  keys_island_subareas_default <- data.frame("subarea"=c("N_SHORE","S_SHORE","N_SHELF","S_SHELF", "WHOLE"))
  ssEnv$keys_sample_groups <-  data.frame("SAMPLE_GROUP"=c("Reference","Control","Case"))

  # filter selected figures, anoamlies and areas passed by user
  figures <- if(is.null(arguments[["figures"]])) keys_figures_default[,1] else arguments$figures
  markers <- if(is.null(arguments[["markers"]]))  keys_markers_default[,1] else arguments$markers

  areas <- if(is.null(arguments[["areas"]]))  keys_areas_default[,1] else arguments$areas

  # check parameters passed by user left areas, figures and anomnalies to work on
  if(sum(figures %in% keys_figures_default[,1])==0)
  {
    message("INFO: ", Sys.time(), " The only allowed figures values are:", keys_figures_default)
    stop("I'm STOPPING HERE!")
  }
  if(sum(markers %in% keys_markers_default[,1])==0)
  {
    message("INFO: ", Sys.time(), " The only allowed markers values are:", keys_markers_default)
    stop("I'm STOPPING HERE!")
  }
  if(sum(areas %in% keys_areas_default[,1])==0)
  {
    message("INFO: ", Sys.time(), " The only allowed areas values are:", keys_areas_default)
    stop("I'm STOPPING HERE!")
  }

  message("INFO: ", Sys.time(), " I will focus on:", paste(markers, collapse = " ", sep =" "), " due to ",  paste(figures, collapse = " ", sep =" "), " of ",  paste(areas, collapse = " ", sep =" "))

  if(!is.null(arguments[["subareas"]]))
  {
    keys_gene_subareas_default <- as.data.frame(keys_gene_subareas_default[keys_gene_subareas_default[,"subarea"]  %in% arguments$subareas,])
    keys_island_subareas_default <- as.data.frame(keys_island_subareas_default[keys_island_subareas_default[,"subarea"] %in% arguments$subareas,])
    if(sum(arguments$subareas %in% c(keys_gene_subareas_default[,1],keys_island_subareas_default))==0)
    {
      message("INFO: ", Sys.time(), " The only allowed areas values are:", keys_areas_default)
      stop("I'm STOPPING HERE!")
    }
  }

  # set working on figures, marker, areas and subareas
  keys_figures <-  data.frame("FIGURE"=figures)
  keys_markers <-  data.frame("MARKER"=markers)
  keys_markers_figures <-  expand.grid("MARKER"=keys_markers[,1],"FIGURE"=keys_figures[,1])

  # create markers figures to work on, cleaning also for BETA only MEAN
  ssEnv$keys_markers_figures <- keys_markers_figures
  levels(ssEnv$keys_markers_figures$FIGURE) <- c(levels(ssEnv$keys_markers_figures$FIGURE),"MEAN")
  ssEnv$keys_markers_figures[ssEnv$keys_markers_figures$MARKER=="BETA","FIGURE"] <- "MEAN"
  ssEnv$keys_markers_figures <- unique(ssEnv$keys_markers_figures)
  ssEnv$keys_markers_figures$COMBINED <- paste0(c(ssEnv$keys_markers_figures$MARKER,ssEnv$keys_markers_figures$FIGURE), collapse ="_")

  keys_areas <- data.frame("AREA"=areas)

  keys_areas_subareas <- data.frame("AREA"="","SUBAREA"="","COMBINED"="")
  keys_areas_subareas <- keys_areas_subareas[-1,]
  ssEnv$keys_areas_subareas_markers_figures <- data.frame("AREA"="","SUBAREA"="","MARKER"="","FIGURE"="")
  ssEnv$keys_areas_subareas_markers_figures <- ssEnv$keys_areas_subareas_markers_figures[-1,]

  if("ISLAND" %in% areas)
  {
    keys_areas_subareas <-  rbind(keys_areas_subareas,expand.grid("AREA"="ISLAND","SUBAREA"=keys_island_subareas_default[,1], "COMBINED"=""))
    ssEnv$keys_areas_subareas_markers_figures <- rbind(ssEnv$keys_areas_subareas_markers_figures,expand.grid("AREA"="ISLAND","SUBAREA"= keys_island_subareas_default[,1],"MARKER"=markers,"FIGURE"=figures))
    message("INFO: ", Sys.time(), " These island's areas will be investigated:", paste(keys_island_subareas_default[,1], collapse = " ", sep=" "))
  }

  if("GENE" %in% areas)
  {
    keys_areas_subareas <- rbind(keys_areas_subareas,expand.grid("AREA"="GENE","SUBAREA"=keys_gene_subareas_default[,1], "COMBINED"=""))
    ssEnv$keys_areas_subareas_markers_figures <- rbind(ssEnv$keys_areas_subareas_markers_figures,expand.grid("AREA"="GENE","SUBAREA"= keys_gene_subareas_default[,1],"MARKER"=markers,"FIGURE"=figures))
    message("INFO: ", Sys.time(), " These gene's areas will be investigated:", paste(keys_gene_subareas_default[,1], collapse = " ", sep=" "))
  }

  if ("DMR" %in% areas)
  {
    keys_areas_subareas <- rbind( data.frame("AREA"="DMR","SUBAREA"="WHOLE", "COMBINED"="DMR_WHOLE"), keys_areas_subareas)
    ssEnv$keys_areas_subareas_markers_figures <- rbind(ssEnv$keys_areas_subareas_markers_figures,expand.grid("AREA"="DMR","SUBAREA"="WHOLE","MARKER"=markers,"FIGURE"=figures))
  }

  if ("CHR" %in% areas)
  {
    keys_areas_subareas <- rbind( data.frame("AREA"="CHR","SUBAREA"="WHOLE", "COMBINED"="CHR_WHOLE"), keys_areas_subareas)
    ssEnv$keys_areas_subareas_markers_figures <- rbind(ssEnv$keys_areas_subareas_markers_figures,expand.grid("AREA"="CHR","SUBAREA"="WHOLE" ,"MARKER"=markers,"FIGURE"=figures))
  }

  if("PROBE" %in% areas)
  {
    keys_areas_subareas <- rbind( data.frame("AREA"="PROBE","SUBAREA"="WHOLE", "COMBINED"="PROBE_WHOLE"), keys_areas_subareas)
    ssEnv$keys_areas_subareas_markers_figures <- rbind(ssEnv$keys_areas_subareas_markers_figures,expand.grid("AREA"="PROBE","SUBAREA"="WHOLE" ,"MARKER"=markers,"FIGURE"=figures))
  }

  combine_not_empty <- function(x)
  {
    paste0(x[x!=""], collapse = "_")
  }

  # force the only FIGURE of BETA as MEAN
  levels(ssEnv$keys_areas_subareas_markers_figures$FIGURE) <- c( levels(ssEnv$keys_areas_subareas_markers_figures$FIGURE),"MEAN")
  ssEnv$keys_areas_subareas_markers_figures$FIGURE[ssEnv$keys_areas_subareas_markers_figures$MARKER=="BETA"] <-"MEAN"
  ssEnv$keys_areas_subareas_markers_figures <- unique(ssEnv$keys_areas_subareas_markers_figures)
  ssEnv$keys_areas_subareas_markers_figures$COMBINED <- apply(ssEnv$keys_areas_subareas_markers_figures[,c("MARKER","FIGURE","AREA","SUBAREA")], 1, combine_not_empty )

  ssEnv$keys_areas_subareas <- keys_areas_subareas
  ssEnv$keys_areas_subareas$COMBINED <- apply(ssEnv$keys_areas_subareas[,c("AREA","SUBAREA")], 1, combine_not_empty )

  # force the only FIGURE of BETA as MEAN
  ssEnv$keys_markers_figures <- keys_markers_figures
  ssEnv$keys_markers_figures <- merge(ssEnv$keys_markers_figures,keys_markers_default, by="MARKER")
  levels(ssEnv$keys_markers_figures$FIGURE) <- c( levels(ssEnv$keys_markers_figures$FIGURE),"MEAN")
  ssEnv$keys_markers_figures[ssEnv$keys_markers_figures$MARKER=="BETA","FIGURE"] <- "MEAN"
  ssEnv$keys_markers_figures$COMBINED <- apply(ssEnv$keys_markers_figures[,c("MARKER","FIGURE")], 1, combine_not_empty )

  ssEnv$keys_areas <- unique(ssEnv$keys_areas_subareas_markers_figures[,"AREA"])

  ssEnv$functionToExport <- c( "analyze_single_sample","deltar_single_sample",
                               "dump_sample_as_bed_file", "delta_single_sample","dir_check_and_create",
                               "file_path_build","analyze_single_sample_both",
                               "sort_by_chr_and_start", "test_match_order", "lesions_get",
                               "mutations_get")


  # to manage progress bar
  if(ssEnv$showprogress)
  {
    handler_settings <- progressr::handlers()
    if ((exists("progress", mode = "function", inherits = TRUE)))
    {
      progressr::handlers(global = TRUE)
      progressr::handlers("progress")
    }
  }

  update_session_info(ssEnv)

  return(ssEnv)
}
