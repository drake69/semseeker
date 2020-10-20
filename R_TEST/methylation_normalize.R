
epimut_normalize_champ <- function(sampleFolder, resultFolder, methodNormalization, useMvalue = FALSE) {
    # if(endsWith(sampleFolder, suffix ='/')) { message('Sample folder must not end with slash!') return() }

    if (.Platform$OS.type == "windows")
        withAutoprint({
            memory.size()
            memory.size(TRUE)
            memory.limit(16000)
        })

    start.time <- Sys.time()

    # IMPORTAZIONE DATI

    # library(parallel)
    library(ChAMP)
    message("Normalizationwarm up ", Sys.time())
    loadMethod <- if (methodNormalization == "SWAN") {
        "minfi"
    } else {
        "ChAMP"
    }
    dataToNormalize <- champ.load(directory = sampleFolder, method = loadMethod)

    normalizationDataFolder <- resultFolder
    if (normalizationDataFolder != "" && !dir.exists(normalizationDataFolder)) {
        dir.create(normalizationDataFolder)
    }

    # browser()

    available.cores <- min(detectCores(all.tests = FALSE, logical = TRUE) - 1, 8)

    # NORMALIZZAZIONE
    message("Starting normalization ", Sys.time())
    normalizedData <- champ.norm(beta = dataToNormalize$beta, rgSet = dataToNormalize$rgSet, mset = dataToNormalize$mset, resultsDir = paste(normalizationDataFolder, "/CHAMP_Normalization_by_", methodNormalization, "/", sep = ""),
        method = methodNormalization, plotBMIQ = FALSE)

    if (any(is.na(normalizedData)))
        stop("Normalized data Produced with NaN values ", Sys.time())

    message("Data normalized ", Sys.time())

    probesNameNormalizedData <- row.names(normalizedData)
    probeFeaturesValue <- data.frame(probe.features)
    probeFeaturesName <- row.names(probeFeaturesValue)

    probeToDiscard <- setdiff(probeFeaturesName, probesNameNormalizedData)
    probeFeaturesValue <- probeFeaturesValue[!row.names(probeFeaturesValue) %in% probeToDiscard, ]
    probeFeaturesValue <- data.frame(Probe = row.names(probeFeaturesValue), probeFeaturesValue)
    probeFeaturesValue$CHR <- droplevels(probeFeaturesValue$CHR)

    probeFeaturesValue <- probeFeaturesValue[order(probeFeaturesValue[, "Probe"], decreasing = FALSE), ]

    names(probeFeaturesValue)[names(probeFeaturesValue) == "MAPINFO"] <- "START"

    normalizedData <- normalizedData[order(row.names(normalizedData), decreasing = FALSE), ]

    if (!test_match_order(row.names(normalizedData), row.names(probeFeaturesValue)))
        stop("Wrong order matching Probes and Mutation!", Sys.time())

    if (!test_match_order(row.names(normalizedData), probeFeaturesValue$Probe))
        stop("Wrong order matching Probes and Mutation!", Sys.time())

    if (useMvalue) {
        normalizedData <- log2(normalizedData/(1 - normalizedData))
    }
    # browser()
    normalizedData <- data.frame(row.names = row.names(probeFeaturesValue), normalizedData)
    write.table(normalizedData, paste(normalizationDataFolder, "/CHAMP_Normalization_by_", methodNormalization, "/", "Normalized_Data.txt", sep = ""), sep = "\t", row.names = TRUE)

    # browser()
    sampleSheet <- if (methodNormalization == "SWAN") {
        data.frame(dataToNormalize$pd@listData)
    } else {
        dataToNormalize$pd
    }

    result <- list(probeFeatures = probeFeaturesValue, sampleSheet = sampleSheet, normalizedData = normalizedData)
    end.time <- Sys.time()
    time.taken <- end.time - start.time
    message("Time to normalize ", time.taken)

    return(result)
}


epimut_normalize_rnbeads <- function(sampleFolder, resultFolder, methodNormalization) {
    # if(endsWith(sampleFolder, suffix ='/')) { message('Sample folder must not end with slash!') return() }

    if (.Platform$OS.type == "windows")
        withAutoprint({
            memory.size()
            memory.size(TRUE)
            memory.limit(16000)
        })

    start.time <- Sys.time()

    # IMPORTAZIONE DATI

    # library(parallel)
    library(ChAMP)
    message("Normalizationwarm up ", Sys.time())
    loadMethod <- if (methodNormalization == "SWAN") {
        "minfi"
    } else {
        "ChAMP"
    }
    dataToNormalize <- champ.load(directory = sampleFolder, method = loadMethod)

    normalizationDataFolder <- sampleFolder
    if (normalizationDataFolder != "" && !dir.exists(normalizationDataFolder)) {
        dir.create(normalizationDataFolder)
    }

    # browser()

    # NORMALIZZAZIONE
    message("Starting normalization ", Sys.time())
    normalizedData <- champ.norm(beta = dataToNormalize$beta, rgSet = dataToNormalize$rgSet, mset = dataToNormalize$mset, resultsDir = paste(normalizationDataFolder, "/CHAMP_Normalization_by_", methodNormalization, "/", sep = ""),
        method = methodNormalization, plotBMIQ = FALSE)

    if (any(is.na(normalizedData)))
        stop("Normalized data Produced with NaN values ", Sys.time())

    message("Data normalized ", Sys.time())

    probesNameNormalizedData <- row.names(normalizedData)
    probeFeaturesValue <- data.frame(probe.features)
    probeFeaturesName <- row.names(probeFeaturesValue)

    probeToDiscard <- setdiff(probeFeaturesName, probesNameNormalizedData)
    probeFeaturesValue <- probeFeaturesValue[!row.names(probeFeaturesValue) %in% probeToDiscard, ]
    probeFeaturesValue <- data.frame(Probe = row.names(probeFeaturesValue), probeFeaturesValue)
    probeFeaturesValue$CHR <- droplevels(probeFeaturesValue$CHR)

    probeFeaturesValue <- probeFeaturesValue[order(probeFeaturesValue[, "Probe"], decreasing = FALSE), ]

    names(probeFeaturesValue)[names(probeFeaturesValue) == "MAPINFO"] <- "START"

    normalizedData <- normalizedData[order(row.names(normalizedData), decreasing = FALSE), ]

    if (!test_match_order(row.names(normalizedData), row.names(probeFeaturesValue)))
        stop("Wrong order matching Probes and Mutation!", Sys.time())

    if (!test_match_order(row.names(normalizedData), probeFeaturesValue$Probe))
        stop("Wrong order matching Probes and Mutation!", Sys.time())

    # browser()
    normalizedData <- data.frame(row.names = row.names(probeFeaturesValue), normalizedData)
    write.table(normalizedData, paste(normalizationDataFolder, "/CHAMP_Normalization_by_", methodNormalization, "/", "Normalized_Data.txt", sep = ""), sep = "\t", row.names = TRUE)


    # browser() NORMALIZZAZIONE
    message("Starting normalization ", Sys.time())
    normalizedData2 <- champ.norm(beta = dataToNormalize$beta, rgSet = dataToNormalize$rgSet, mset = dataToNormalize$mset, resultsDir = paste(normalizationDataFolder, "/CHAMP_Normalization_by_", methodNormalization, "/", sep = ""),
        method = methodNormalization, plotBMIQ = FALSE)

    if (any(is.na(normalizedData2)))
        stop("Normalized data Produced with NaN values ", Sys.time())

    message("Data normalized ", Sys.time())

    probesNamenormalizedData2 <- row.names(normalizedData2)
    probeFeaturesValue <- data.frame(probe.features)
    probeFeaturesName <- row.names(probeFeaturesValue)

    probeToDiscard <- setdiff(probeFeaturesName, probesNamenormalizedData2)
    probeFeaturesValue <- probeFeaturesValue[!row.names(probeFeaturesValue) %in% probeToDiscard, ]
    probeFeaturesValue <- data.frame(Probe = row.names(probeFeaturesValue), probeFeaturesValue)
    probeFeaturesValue$CHR <- droplevels(probeFeaturesValue$CHR)

    probeFeaturesValue <- probeFeaturesValue[order(probeFeaturesValue[, "Probe"], decreasing = FALSE), ]

    names(probeFeaturesValue)[names(probeFeaturesValue) == "MAPINFO"] <- "START"

    normalizedData2 <- normalizedData2[order(row.names(normalizedData2), decreasing = FALSE), ]

    if (!test_match_order(row.names(normalizedData2), row.names(probeFeaturesValue)))
        stop("Wrong order matching Probes and Mutation!", Sys.time())

    if (!test_match_order(row.names(normalizedData2), probeFeaturesValue$Probe))
        stop("Wrong order matching Probes and Mutation!", Sys.time())

    # browser()
    normalizedData2 <- data.frame(row.names = row.names(probeFeaturesValue), normalizedData2)
    write.table(normalizedData2, paste(normalizationDataFolder, "/CHAMP_Normalization_by_", methodNormalization, "/", "Normalized_Data_2.txt", sep = ""), sep = "\t", row.names = TRUE)

    check <- normalizedData != normalizedData2
    print(table(check))

    browser()
    sampleSheet <- if (methodNormalization == "SWAN") {
        data.frame(dataToNormalize$pd@listData)
    } else {
        dataToNormalize$pd
    }

    result <- list(probeFeatures = probeFeaturesValue, sampleSheet = sampleSheet, normalizedData = normalizedData)
    end.time <- Sys.time()
    time.taken <- end.time - start.time
    message("Time to normalize ", time.taken)

    return(result)
}
