test_that("test_init_env", {

  # run this test as first other wise the reuse of session is not testable
  tempFolder <- tempFolders[1]
  tempFolders <- tempFolders[-1]
  ssEnv <- semseeker:::init_env(tempFolder, parallel_strategy = parallel_strategy, areas = c("GENE"), subareas= c("BODY"), showprogress=TRUE, start_fresh=TRUE, maxResources = 10)
  semseeker:::close_env()

  ssEnv <- semseeker:::init_env(tempFolder, parallel_strategy = parallel_strategy,areas = c("GENE","DMR"), showprogress=TRUE, start_fresh=TRUE)
  testthat::expect_true(length(unique(ssEnv$keys_areas_subareas[,"AREA"]))!=1)

  unlink(tempFolder,recursive = TRUE)
  assign("ssEnv", NULL, envir=.pkgglobalenv)

  ####################################################################################

  testthat::expect_error( semseeker:::init_env(tempFolder, parallel_strategy = parallel_strategy, figures="HYPPO"), "ERROR:")
  ####################################################################################


  testthat::expect_error( semseeker:::init_env(tempFolder, parallel_strategy = parallel_strategy, markers="HYPPO"), "ERROR:")

  ####################################################################################

  testthat::expect_error( semseeker:::init_env(tempFolder, parallel_strategy = parallel_strategy, areas="HYPPO"), "ERROR:")

  ####################################################################################

  # expect_error( semseeker:::init_env(tempFolder, parallel_strategy ="cluster"), "ERROR:")

  ####################################################################################

  testthat::expect_error( semseeker:::init_env(tempFolder, parallel_strategy = parallel_strategy, subareas="HYPPO"), "ERROR:")


  testthat::expect_error( semseeker:::init_env(tempFolder, parallel_strategy = parallel_strategy, test_it="HYPPO"), "ERROR:")

  ####################################################################################

  unlink(tempFolder,recursive = TRUE)
  # semseeker:::close_env()
}
)
