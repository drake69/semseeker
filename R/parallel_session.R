parallel_session <- function()
{
  #
  ssEnv <- get_session_info()
  parallel_strategy <- ssEnv$parallel_strategy

  # check if os is macos
  if (Sys.info()["sysname"] == "Darwin" & parallel_strategy != "sequential")
  {
    env_var <- Sys.getenv("OBJC_DISABLE_INITIALIZE_FORK_SAFETY")
    if (env_var != "YES")
    {
      log_event("ERROR: Setting OBJC_DISABLE_INITIALIZE_FORK_SAFETY must be YES to work in multiprocess. \n
        execute: export OBJC_DISABLE_INITIALIZE_FORK_SAFETY=YES at shell or
         Sys.setenv(OBJC_DISABLE_INITIALIZE_FORK_SAFETY='YES') ! ")
      # Sys.setenv(OBJC_DISABLE_INITIALIZE_FORK_SAFETY="YES")
      stop()
    }
  }

  if(parallelly::supportsMulticore())
    options(parallelly.fork.enable= TRUE)
  else
    options(parallelly.fork.enable= FALSE)

  chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")

  if (nzchar(chk) && chk == "TRUE") {
    # use 2 cores in CRAN/Travis/AppVeyor
    nCore <- 2L
  } else {
    # use all cores in devtools::test()
    nCore <- future::availableCores() - 1
    nCore <- if(floor(future::availableCores() * ssEnv$maxResources/100 ) > nCore ) nCore else floor(future::availableCores() * ssEnv$maxResources/100 )
  }
  # permutation cluster
  outFile <- file.path(ssEnv$session_folder, "cluster_r.out")

  #
  ssEnv$parallel <- data.frame("parallel_strategy"="")
  ssEnv$parallel$parallel_strategy <- parallel_strategy
  ssEnv$parallel$nCore <- nCore

  options(doFuture.foreach.export = ".export-and-automatic-with-warning")
  # doFuture
  # is changed in onload see zzz.R
  # options(parallelly.fork.enable= FALSE)
  # options(future.globals.resolve = TRUE)
  #allow export of object of 32gb with future
  options(future.globals.maxSize= 32 * 1024^3)
  doFuture::registerDoFuture()

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
    if (!is.null(ssEnv$cluster_workers))
    {
      future::plan( future::cluster, workers = ssEnv$cluster_workers)
      log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), " I will work with a cluster with:",ssEnv$cluster_workers)
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
    options(parallelly.fork.enable= FALSE)
    future::plan(strategy = future::sequential)
    log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), " I will work in sequential mode")
  }

  update_session_info(ssEnv)

}
