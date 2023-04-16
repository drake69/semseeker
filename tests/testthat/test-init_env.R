test_that("test_init_env", {


  tmp <- tempdir()
  tempFolder <- paste(tmp,"/semseeker/",stringi::stri_rand_strings(1, 7, pattern = "[A-Za-z0-9]"),sep="")

  ####################################################################################

  expect_error(init_env(tempFolder, parallel_strategy = parallel_strategy, figures="HYPPO"), "I'm STOPPING HERE!")

  ####################################################################################

  expect_error(init_env(tempFolder, parallel_strategy = parallel_strategy, anomalies="HYPPO"), "I'm STOPPING HERE!")

  ####################################################################################

  expect_error(init_env(tempFolder, parallel_strategy = parallel_strategy, metaareas="HYPPO"), "I'm STOPPING HERE!")

  ####################################################################################

  expect_error(init_env(tempFolder, parallel_strategy ="cluster"), "I'm STOPPING HERE!")

  ####################################################################################

  expect_error(init_env(tempFolder, parallel_strategy = parallel_strategy, metaareassub="HYPPO"), "I'm STOPPING HERE!")

  ####################################################################################
  # ssEnv <- init_env(tempFolder, parallel_strategy = parallel_strategy, showprogress=TRUE)
  # testthat::expect_true(ssEnv$shoprogress==TRUE)

  unlink(tempFolder,recursive = TRUE)
  # close_env()

}
)
