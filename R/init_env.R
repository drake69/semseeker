#' init ssEnvonment
#'
#' @param result_folder where result of semseeker will bestored
#' @param parallel_strategy which strategy to use for parallel executio see future vignete: possibile values, none, multisession,sequential, multicore, cluster
#' @param maxResources percentage of how many available cores will be used default 90 percent, rounded to the lowest integer
#' @param ... other options to filter elaborations
#'
#' @return the working ssEnvonment
init_env <- function(result_folder, maxResources = 90, ...)
{
  tryCatch(
    {
      test_it <- list(...)
    },
    error = function(cond)  {
      log_event ("ERROR: ", format(Sys.time(), "%a %b %d %X %Y"), " Function's arguments must be passed explicitily !")
      log_event(cond)
      stop()
    }
  )

  arguments <- list(...)

  # remove all spaces from all items of arguments
  arguments <- lapply(arguments, function(x) gsub(" ", "", x))

  start_fresh <- TRUE
  if(!is.null(arguments[["start_fresh"]]))
    start_fresh <- arguments$start_fresh
  arguments[["start_fresh"]] <- NULL

  if(start_fresh)
    unlink(result_folder, recursive = TRUE, force = TRUE)

  #allow export of object of 32gb with future
  options(future.globals.maxSize= 32 * 1024^3)

  ssEnv <- get_session_info(result_folder)

  # utils::data("PROBES")
  # utils::data("PROBES_CHR_CHR")
  ssEnv$seed <- 7658776
  set.seed(7658776)
  ssEnv$session_folder <-  dir_check_and_create(result_folder,c("Log"))

  log_event("INFO: ##############################################################################")
  log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), " Job Started !")

  arguments <- set_env_variable(ssEnv, arguments,"alpha",0.05)
  ssEnv <- get_session_info(result_folder)
  arguments <- set_env_variable(ssEnv, arguments,"verbosity",1)
  ssEnv <- get_session_info(result_folder)
  arguments <- set_env_variable(ssEnv, arguments,"sex_chromosome_remove",FALSE)
  ssEnv <- get_session_info(result_folder)
  arguments <- set_env_variable(ssEnv, arguments,"opencl",FALSE)
  ssEnv <- get_session_info(result_folder)
  arguments <- set_env_variable(ssEnv, arguments,"bonferroni_threshold",0.05)
  ssEnv <- get_session_info(result_folder)
  arguments <- set_env_variable(ssEnv, arguments,"iqrTimes",3)
  ssEnv <- get_session_info(result_folder)
  arguments <- set_env_variable(ssEnv, arguments,"sliding_window_size",11)
  ssEnv <- get_session_info(result_folder)
  arguments <- set_env_variable(ssEnv, arguments,"epiquantile",4)
  ssEnv <- get_session_info(result_folder)
  arguments <- set_env_variable(ssEnv, arguments,"maxResources",90)
  ssEnv <- get_session_info(result_folder)
  arguments <- set_env_variable(ssEnv, arguments,"parallel_strategy","sequential")
  ssEnv <- get_session_info(result_folder)
  arguments <- set_env_variable(ssEnv, arguments,"tech","")
  ssEnv <- get_session_info(result_folder)
  arguments <- set_env_variable(ssEnv, arguments,"showprogress",FALSE)
  ssEnv <- get_session_info(result_folder)
  arguments <- set_env_variable(ssEnv, arguments,"signal_intrasample",FALSE)
  ssEnv <- get_session_info(result_folder)
  original_colors <- c('#b9e192', '#b3c7f7', '#f8b8d0','#f194b8','#b9ef92', '#ffefb6', '#cfebb6')
  arguments <- set_env_variable(ssEnv, arguments,"color_palette",original_colors)
  ssEnv <- get_session_info(result_folder)
  darker_colors <- grDevices::adjustcolor(original_colors, alpha.f = 0.5)
  darker_colors <- c("blue","red","purple","green","yellow","orange","brown")
  arguments <- set_env_variable(ssEnv, arguments,"color_palette_darker",darker_colors)
  ssEnv <- get_session_info(result_folder)
  arguments <- set_env_variable(ssEnv, arguments,"cluster_workers",NULL)
  ssEnv <- get_session_info(result_folder)

  dry_run <- FALSE
  if(!is.null(arguments[["dry_run"]]))
    dry_run <- arguments$dry_run
  arguments[["dry_run"]] <- NULL
  if(dry_run)
    ssEnv$verbosity <- 4

  tmp <- tempdir()
  log_event("INFO: ",format(Sys.time(), "%a %b %d %X %Y")," data will saved in this folder:", result_folder)
  ssEnv$temp_folder <-  paste(tmp,"/semseeker/",stringi::stri_rand_strings(1, 7, pattern = "[A-Za-z0-9]"),sep="")
  ssEnv$result_folderData <-  dir_check_and_create(result_folder, "Data")
  ssEnv$result_folderChart <-    dir_check_and_create(result_folder, "Chart")
  ssEnv$result_folderInference <-    dir_check_and_create(result_folder, "Inference")
  ssEnv$result_folderPathway <-    dir_check_and_create(result_folder, "Pathway")
  ssEnv$result_folderPhenotype <-    dir_check_and_create(result_folder, "Phenotype")
  ssEnv$result_folderEuristic <-  dir_check_and_create(result_folder,"Euristic")
  ssEnv$session_folder <-  dir_check_and_create(result_folder,c("Log"))
  random_file_name <- paste(stringi::stri_rand_strings(1, 7, pattern = "[A-Za-z0-9]"),".log", sep="")

  if (sink.number() != 0)
    sink(NULL)
  sink(file.path(ssEnv$session_folder,"session_output.log"), split = TRUE, append = TRUE)

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
  # permutation cluster
  outFile <- file.path(ssEnv$session_folder, "cluster_r.out")

  options(doFuture.foreach.export = ".export-and-automatic-with-warning")
  doFuture::registerDoFuture()


  parallel_strategy <- ssEnv$parallel_strategy
  # TODO: improve planning parallel management using also cluster
  if(parallel_strategy=="multisession")
  {
    future::plan( future::multisession, workers = nCore)
    log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), " I will work in multisession with:", nCore, " Cores")
  }
  if(parallel_strategy=="multicore")
  {
    future::plan( future::multicore, workers = nCore)
    log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), " I will work in muticore with:", nCore," Cores")
  }
  if(parallel_strategy=="cluster")
  {
    # log_event ("ERROR: ", format(Sys.time(), "%a %b %d %X %Y"), " Cluster feature not implemented!")
    # stop("I'm STOPPING HERE!")
    if (!is.null(ssEnv$cluster_workers))
    {
      future::plan( future::cluster, workers = ssEnv$cluster_workers)
      log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), " I will work with a cluster with:",ssEnv$cluster_workers)
      # ssEnv$showprogress <- FALSE
    }
    else
    {
      future::plan( future::cluster, workers = nCore)
      log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), " I will work with a cluster with:", nCore," Cores")
    }
  }
  if(parallel_strategy!="multisession" & parallel_strategy!="multicore"
     & parallel_strategy!="cluster")
  {
    future::plan( future::sequential)
    log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), " I will work in sequential mode")
  }


  # ssEnvTemp <- get_session_info(result_folder)
  # if (!is.null(ssEnvTemp) & !length(ssEnvTemp)<2)
  # {
  #   # we had to return old session info this init could be called by other than semseeker function, it means
  #   # we don't know which marker or figure were identified
  #   log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), " Reusing old session info !")
  #   log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), " I will focus on: ", paste(unique(ssEnvTemp$keys_markers_figures[,"MARKER"]), collapse = " ", sep =" "), " due to ",
  #     paste(unique(ssEnvTemp$keys_markers_figures[,"FIGURE"]), collapse = " ", sep =" "),
  #     " of ",  paste(unique(ssEnvTemp$keys_areas_subareas[,"AREA"]), collapse = " ", sep =" "),
  #     " for ",  paste(unique(ssEnvTemp$keys_areas_subareas[,"SUBAREA"]), collapse = " ", sep =" "))
  #   if(dir.exists(ssEnv$session_folder))
  #   {
  #     update_session_info(ssEnvTemp)
  #     return(ssEnvTemp)
  #   }
  # }

  # set default values
  keys_figures_default <-  data.frame("FIGURE"=c("HYPO", "HYPER", "BOTH","BOTHSUM"))

  ssEnv$keys_figures_default <- keys_figures_default

  keys_markers_default <-  data.frame("MARKER"=c("MUTATIONS","LESIONS","DELTAS","DELTAQ","DELTAR","DELTARQ","SIGNAL"))
  keys_markers_default$SUFFIX <-  c("","","","","","","")
  # keys_markers_default$SUFFIX <-  c("","","","", ssEnv$epiquantile ,"",ssEnv$epiquantile)
  keys_markers_default$EXT <-  c("bed","bed","bedgraph","bed","bedgraph","bed","bedgraph")
  ssEnv$keys_markers_default <- keys_markers_default


  keys_markers_figure_default <-  reshape::expand.grid.df(keys_markers_default, keys_figures_default)
  keys_markers_figure_default[keys_markers_figure_default$MARKER=="SIGNAL","FIGURE"] <- "MEAN"
  keys_markers_figure_default <- unique(keys_markers_figure_default)
  ssEnv$keys_markers_figure_default <- keys_markers_figure_default


  keys_areas_default <- data.frame("areas"=c("GENE","ISLAND","DMR","CHR","PROBE"))
  keys_gene_subareas_default <- data.frame("subarea"=c("BODY","TSS1500","TSS200","1STEXON","3UTR","5UTR","EXONBND","WHOLE"))
  keys_island_subareas_default <- data.frame("subarea"=c("N_SHORE","S_SHORE","N_SHELF","S_SHELF", "WHOLE"))

  ssEnv$keys_sample_groups <-  data.frame("SAMPLE_GROUP"=c("Reference","Control","Case"))

  # filter selected figures, anoamlies and areas passed by user
  figures <- if(is.null(arguments[["figures"]])) keys_figures_default[,1] else arguments$figures
  arguments[["figures"]] <- NULL

  markers <- if(is.null(arguments[["markers"]]))  keys_markers_default[,1] else arguments$markers
  arguments[["markers"]] <- NULL

  areas <- if(is.null(arguments[["areas"]]))  keys_areas_default[,1] else unique(arguments$areas)
  arguments[["areas"]] <- NULL

  # check parameters passed by user left areas, figures and anomnalies to work on
  if(sum(figures %in% keys_figures_default[,1])==0)
  {
    log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), " The only allowed figures values are:", keys_figures_default)
    stop("I'm STOPPING HERE!")
  }
  if(sum(markers %in% keys_markers_default[,1])==0)
  {
    log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), " The only allowed markers values are:", keys_markers_default)
    stop("I'm STOPPING HERE!")
  }
  if(sum(areas %in% keys_areas_default[,1])==0)
  {
    log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), " The only allowed areas values are:", keys_areas_default)
    stop("I'm STOPPING HERE!")
  }


  if(!is.null(arguments[["subareas"]]))
  {
    keys_gene_subareas_default <- as.data.frame(keys_gene_subareas_default[keys_gene_subareas_default[,"subarea"]  %in% arguments$subareas,])
    keys_island_subareas_default <- as.data.frame(keys_island_subareas_default[keys_island_subareas_default[,"subarea"] %in% arguments$subareas,])
    if(sum(arguments$subareas %in% c(keys_gene_subareas_default[,1],keys_island_subareas_default))==0)
    {
      log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), " The only allowed areas values are:", keys_areas_default)
      stop("I'm STOPPING HERE!")
    }
    arguments[["subareas"]] <- NULL
  }

  # set working on figures, marker, areas and subareas
  keys_figures <-  data.frame("FIGURE"=figures)
  keys_markers <-  data.frame("MARKER"=markers)
  keys_markers_figures <-  expand.grid("MARKER"=keys_markers[,1],"FIGURE"=keys_figures[,1])

  # create markers figures to work on, cleaning also for SIGNAL only MEAN
  ssEnv$keys_markers_figures <- keys_markers_figures
  levels(ssEnv$keys_markers_figures$FIGURE) <- c(levels(ssEnv$keys_markers_figures$FIGURE),"MEAN")
  ssEnv$keys_markers_figures[ssEnv$keys_markers_figures$MARKER=="SIGNAL","FIGURE"] <- "MEAN"
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
    log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), " These island's areas will be investigated:", paste(keys_island_subareas_default[,1], collapse = " ", sep=" "))
  }

  if("GENE" %in% areas)
  {
    keys_areas_subareas <- rbind(keys_areas_subareas,expand.grid("AREA"="GENE","SUBAREA"=keys_gene_subareas_default[,1], "COMBINED"=""))
    ssEnv$keys_areas_subareas_markers_figures <- rbind(ssEnv$keys_areas_subareas_markers_figures,expand.grid("AREA"="GENE","SUBAREA"= keys_gene_subareas_default[,1],"MARKER"=markers,"FIGURE"=figures))
    log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), " These gene's areas will be investigated:", paste(keys_gene_subareas_default[,1], collapse = " ", sep=" "))
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
    # paste0(x[x!=""], collapse = "_")
    # remove spaces
    gsub(" ","",paste0(x[x!=""], collapse = "_"))
  }

  # force the only FIGURE of SIGNAL as MEAN
  levels(ssEnv$keys_areas_subareas_markers_figures$FIGURE) <- c( levels(ssEnv$keys_areas_subareas_markers_figures$FIGURE),"MEAN")
  ssEnv$keys_areas_subareas_markers_figures$FIGURE[ssEnv$keys_areas_subareas_markers_figures$MARKER=="SIGNAL"] <-"MEAN"
  ssEnv$keys_areas_subareas_markers_figures <- unique(ssEnv$keys_areas_subareas_markers_figures)
  ssEnv$keys_areas_subareas_markers_figures$COMBINED <- apply(ssEnv$keys_areas_subareas_markers_figures[,c("MARKER","FIGURE","AREA","SUBAREA")], 1, combine_not_empty )

  ssEnv$keys_areas_subareas <- keys_areas_subareas
  ssEnv$keys_areas_subareas$COMBINED <- apply(ssEnv$keys_areas_subareas[,c("AREA","SUBAREA")], 1, combine_not_empty )

  # force the only FIGURE of SIGNAL as MEAN
  ssEnv$keys_markers_figures <- keys_markers_figures
  ssEnv$keys_markers_figures <- merge(ssEnv$keys_markers_figures,keys_markers_default, by="MARKER")
  levels(ssEnv$keys_markers_figures$FIGURE) <- c( levels(ssEnv$keys_markers_figures$FIGURE),"MEAN")
  ssEnv$keys_markers_figures[ssEnv$keys_markers_figures$MARKER=="SIGNAL","FIGURE"] <- "MEAN"
  ssEnv$keys_markers_figures$COMBINED <- apply(ssEnv$keys_markers_figures[,c("MARKER","FIGURE")], 1, combine_not_empty )

  ssEnv$keys_areas <- unique(ssEnv$keys_areas_subareas_markers_figures[,"AREA"])


  keys <- unique(ssEnv$keys_areas_subareas_markers_figures)
  levels(keys$FIGURE) <- c(levels(keys$FIGURE),"HYPER_HYPO")
  levels(keys$SUBAREA) <- c(levels(keys$SUBAREA),"GENE_PARTS")
  keys[keys$FIGURE!="BOTH" & keys$FIGURE!="MEAN","FIGURE"] <- "HYPER_HYPO"
  keys[keys$SUBAREA!="WHOLE","SUBAREA"] <- "GENE_PARTS"
  keys$COMBINED <- paste(keys$AREA,keys$SUBAREA,keys$MARKER,keys$FIGURE, sep="_")
  keys <- keys[keys$AREA=="GENE",]
  keys <- keys[!duplicated(keys),]
  ssEnv$keys_for_pathway <- keys

  ssEnv$parallel <- data.frame("strategy"=parallel_strategy, "nCore"=nCore)

  ssEnv$functionToExport <- c( "analyze_single_sample","deltar_single_sample",
                               "dump_sample_as_bed_file", "delta_single_sample","dir_check_and_create",
                               "file_path_build","analyze_single_sample_both",
                               "sort_by_chr_and_start", "test_match_order", "lesions_get",
                               "mutations_get")


  # to manage progress bar
  if(ssEnv$showprogress)
  {
    handler_settings <- progressr::handlers()
    # if (!(exists("progress", mode = "function", inherits = TRUE)))
    # {
    #   progressr::handlers(global = TRUE)
    #   progressr::handlers("progress")
    # }
    if (!(exists("cli", mode = "function", inherits = TRUE)))
    {
      # check if handler is already registered
      if(!("cli" %in% handler_settings$handler))
      {
        # check if handlers is on the stack
        if(!("cli" %in% handler_settings$stack))
        {
          progressr::handlers(global = TRUE)
          progressr::handlers("cli")
        }
      }
    }
  }

  update_session_info(ssEnv)
  log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), " I will focus on:", paste(unique(keys_markers_figures$MARKER), collapse = " ", sep =" "), " due to ",  paste(unique(keys_markers_figures$FIGURE), collapse = " ", sep =" "), " of ",  paste(areas, collapse = " ", sep =" "))

  # check length of arguments
  if(length(arguments)!=0)
  {
    log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), " This options are not recognized: ", paste(arguments, collapse = " ", sep =" "))
    stop()
    }


  if(dry_run)
    {
      # out at console as pretty table

      knitr::kable(as.data.frame(ssEnv$keys_areas_subareas_markers_figures), format = "pipe", caption = "Selection:")
      message(ssEnv$keys_areas_subareas_markers_figures)
      stop("INFO: Dry run is requested. Exiting now.")
    }

  return(ssEnv)
}
