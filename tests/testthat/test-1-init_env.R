test_that("test_init_env", {


  tmp <- tempdir()
  tempFolder <- paste(tmp,"/semseeker/",stringi::stri_rand_strings(1, 7, pattern = "[A-Za-z0-9]"),sep="")

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
  # ssEnv <- semseeker:::init_env(tempFolder, parallel_strategy = parallel_strategy, showprogress=TRUE)
  # testthat::expect_true(ssEnv$shoprogress==TRUE)

  unlink(tempFolder,recursive = TRUE)
  # semseeker:::close_env()

}
)
