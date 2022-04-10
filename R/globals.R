# if(getRversion() >= "2.15.1")  utils::globalVariables(c( "ssEnv$resultFolderData","chartFolder","ssEnv$logFolder","computationCluster",
#                                                          "ssEnv$resultFolderChart","ssEnv$resultFolderInference","ssEnv$resultFolderEuristic",
#                                                          "ssEnv$keys_anomalies","ssEnv$keys_figures","ssEnv$keys_populations","ssEnv$keys",
#                                                          "probes","functionToExport"), add = FALSE)

# my_env <- new.env(parent = emptyenv())
# functionToExport <<- ssEnv$resultFolderData <<- chartFolder <<- ssEnv$logFolder <<- computationCluster <<- NULL
# ssEnv$resultFolderChart <<- ssEnv$resultFolderInference <<- ssEnv$resultFolderEuristic <<-NULL
# ssEnv$keys_anomalies <<- ssEnv$keys_figures <<- ssEnv$keys_populations <<- ssEnv$keys <<- probes <<- foreachIndex <<- NULL

  # utils::globalVariables( names= c( "ssEnv$resultFolderData","chartFolder","ssEnv$logFolder","computationCluster",
  #                           "ssEnv$resultFolderChart","ssEnv$resultFolderInference","ssEnv$resultFolderEuristic",
  #                           "ssEnv$keys_anomalies","ssEnv$keys_figures","ssEnv$keys_populations","ssEnv$keys",
  #                           "probes","foreachIndex"), add = FALSE, package = "semseeker")

  # ssEnv <- new.env(parent = emptyenv())
  # ssEnv$resultFolderData <- NULL
  # chartFolder  <- NULL
  # ssEnv$logFolder   <- NULL
  # computationCluster <- NULL
  # ssEnv$resultFolderChart  <- NULL
  # ssEnv$resultFolderInference  <- NULL
  # ssEnv$resultFolderEuristic   <- NULL
  # ssEnv$keys_anomalies  <- NULL
  # ssEnv$keys_figures  <- NULL
  # ssEnv$keys_populations   <- NULL
  # ssEnv$keys   <- NULL
  # probes  <- NULL
  # foreachIndex  <- NULL

# ssEnv <- new.env(parent = emptyenv())

