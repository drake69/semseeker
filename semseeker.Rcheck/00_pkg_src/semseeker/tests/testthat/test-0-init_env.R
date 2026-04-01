test_that("test_init_env", {

  # run this test as first other wise the reuse of session is not testable
  tempFolder <- tempFolders[1]
  tempFolders <- tempFolders[-1]
  ssEnv <- semseeker:::init_env(tempFolder, parallel_strategy = parallel_strategy, areas = c("GENE"), subareas= c("BODY"), showprogress=TRUE)
  ssEnv <- semseeker:::init_env(tempFolder, parallel_strategy = parallel_strategy,areas = c("GENE","DMR"), showprogress=TRUE)
  testthat::expect_true(length(unique(ssEnv$keys_areas_subareas[,"AREA"]))==1)

  unlink(tempFolder,recursive = TRUE)
  assign("ssEnv", NULL, envir=.pkgglobalenv)

  ####################################################################################

  expect_error( semseeker:::init_env(tempFolder, parallel_strategy = parallel_strategy, figures="HYPPO"), "I'm STOPPING HERE!")
  ####################################################################################


  expect_error( semseeker:::init_env(tempFolder, parallel_strategy = parallel_strategy, markers="HYPPO"), "I'm STOPPING HERE!")

  ####################################################################################

  expect_error( semseeker:::init_env(tempFolder, parallel_strategy = parallel_strategy, areas="HYPPO"), "I'm STOPPING HERE!")

  ####################################################################################

  expect_error( semseeker:::init_env(tempFolder, parallel_strategy ="cluster"), "I'm STOPPING HERE!")

  ####################################################################################

  expect_error( semseeker:::init_env(tempFolder, parallel_strategy = parallel_strategy, areas="HYPPO"), "I'm STOPPING HERE!")

  ####################################################################################

  unlink(tempFolder,recursive = TRUE)
  # semseeker:::close_env()
}
)
