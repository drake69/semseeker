# if(getRversion() >= "2.15.1")  utils::globalVariables(c( "resultFolderData","chartFolder","logFolder","computationCluster",
#                                                          "resultFolderChart","resultFolderInference","resultFolderEuristic",
#                                                          "keys_anomalies","keys_figures","keys_populations","keys",
#                                                          "probes"), add = FALSE)

# my_env <- new.env(parent = emptyenv())

resultFolderData <-chartFolder <- logFolder <- computationCluster <- NULL
resultFolderChart <- resultFolderInference<- resultFolderEuristic <-NULL
keys_anomalies <- keys_figures <- keys_populations<- keys <- probes <- foreachIndex <- NULL
