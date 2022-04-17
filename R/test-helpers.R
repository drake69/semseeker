# if(getRversion() >= "2.15.1")  utils::globalVariables(c( "ssEnv$result_folderData","chartFolder","ssEnv$logFolder","computationCluster",
#                                                          "ssEnv$result_folderChart","ssEnv$result_folderInference","ssEnv$result_folderEuristic",
#                                                          "ssEnv$keys_anomalies","ssEnv$keys_figures","ssEnv$keys_populations","ssEnv$keys",
#                                                          "probes"), add = FALSE)

# my_env <- new.env(parent = emptyenv())
# functionToExport <<- ssEnv$result_folderData <<- chartFolder <<- ssEnv$logFolder <<- computationCluster <<- NULL
# ssEnv$result_folderChart <<- ssEnv$result_folderInference <<- ssEnv$result_folderEuristic <<-NULL
# ssEnv$keys_anomalies <<- ssEnv$keys_figures <<- ssEnv$keys_populations <<- ssEnv$keys <<- probes <<- foreachIndex <<- NULL

# utils::globalVariables( names= c( "ssEnv$result_folderData","chartFolder","ssEnv$logFolder","computationCluster",
#                           "ssEnv$result_folderChart","ssEnv$result_folderInference","ssEnv$result_folderEuristic",
#                           "ssEnv$keys_anomalies","ssEnv$keys_figures","ssEnv$keys_populations","ssEnv$keys",
#                           "probes","foreachIndex"), add = FALSE, package = "semseeker")

 # ssEnv <- new.env(parent = emptyenv())
# ssEnv$result_folderData <- NULL
# chartFolder  <- NULL
# ssEnv$logFolder   <- NULL
# computationCluster <- NULL
# ssEnv$result_folderChart  <- NULL
# ssEnv$result_folderInference  <- NULL
# ssEnv$result_folderEuristic   <- NULL
# ssEnv$keys_anomalies  <- NULL
# ssEnv$keys_figures  <- NULL
# ssEnv$keys_populations   <- NULL
# ssEnv$keys   <- NULL
# probes  <- NULL
# foreachIndex  <- NULL
#


