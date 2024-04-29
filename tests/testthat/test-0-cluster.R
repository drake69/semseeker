test_that("cluster", {

  # source("./tests/testthat/setup.R")
  tempFolder <- "~/Documenti/tmp/test_1"
  ssEnv <- semseeker:::init_env(tempFolder, parallel_strategy = parallel_strategy,
    areas = c("GENE"), subareas= c("BODY"), showprogress=FALSE, start_fresh=FALSE,
    cluster_workers="10.161.5.13")
  library(future)
  # gg %<-% {
  #   # p <- ggplot(data = d, aes(x = carat, y = price)) +
  #   #   geom_point(aes(text = paste("Clarity:", clarity)), size = 4) +
  #   #   geom_smooth(aes(colour = cut, fill = cut)) + facet_wrap(~ cut)
  #   # ggplotly(p)
  #   semseeker:::get_session_info("/home/lcorsaro/Desktop/test")
  # }
  #
  # gg

  gg1 %<-% {
    # p <- ggplot(data = d, aes(x = carat, y = price)) +
    #   geom_point(aes(text = paste("Clarity:", clarity)), size = 4) +
    #   geom_smooth(aes(colour = cut, fill = cut)) + facet_wrap(~ cut)
    # ggplotly(p)
    semseeker:::get_session_info()
  }
  testthat::expect_true(gg1$cluster_workers=="10.161.5.13")
  testthat::expect_true(gg1$result_folderData==paste(tempFolder,"/Data",sep=""))
})
