#' for each sample given the mutations return the lesions
#'
#' @param mutationAnnotatedSorted dataframe of mutations annotated and sorted
#' per the input sample
#' @param slidingWindowSize window size to calculate the hypergeometric
#' distribution
#' @param sampleName name of the sample which working with
#' @param bonferroniThreshold accepted value to consider pValue valid
#'
#' @return data frame with lesioned/not lesioned probes identified by 1/0
#' @export
#'
getLesions <- function(mutationAnnotatedSorted, slidingWindowSize, sampleName, bonferroniThreshold = 0.05, probeFeatures) {
  LESIONS <- NULL

  colnames_required <- c("CHR", "START", "END")
  checkr::check_colnames(probeFeatures, colnames_required, error = TRUE)
  checkr::check_colnames(mutationAnnotatedSorted, "MUTATIONS", error = TRUE)

  browser()
  mutationAnnotatedSortedLocal <- sortByCHRandSTART(mutationAnnotatedSorted)

  if (!test_match_order(row.names(mutationAnnotatedSortedLocal), mutationAnnotatedSortedLocal$PROBE)) {
    stop("Wrong order matching Probes and Mutations annotated!", Sys.time())
  }

  mutationAnnotatedSortedLocal$CHR <- droplevels(mutationAnnotatedSortedLocal$CHR)

  if (length(mutationAnnotatedSortedLocal$CHR) <= slidingWindowSize) {
    browser()
  }

  missedWindowLength <- ((slidingWindowSize - 1) / 2)

  # browser()
  startingValue <- mutationAnnotatedSortedLocal[1:missedWindowLength,]
  startingValue$MUTATIONS <- rep(0,missedWindowLength)

  endingValue <- mutationAnnotatedSortedLocal[(length(mutationAnnotatedSortedLocal) - missedWindowLength):length(mutationAnnotatedSortedLocal),]
  endingValue$MUTATIONS <- rep(0,missedWindowLength)

  # browser()
  mutationAnnotatedSortedLocal <- rbind(startingValue, mutationAnnotatedSortedLocal, endingValue)

  mutationAnnotatedSortedWindowed <- zoo::zoo(x = mutationAnnotatedSortedLocal[, "MUTATIONS"])
  message(sampleName, " ", "Got mutationAnnotatedSortedWindowed ", Sys.time())

  sumIntoCentreWindow <- function(x) {
    zoo::rollapply(x, width = slidingWindowSize, by = 1, FUN = sum, align = "center")
  }

  ## check burden into the window (sum)
  mutationAnnotatedSortedWindowedSum <- sumIntoCentreWindow(mutationAnnotatedSortedWindowed$MUTATIONS)

  message(sampleName, " ", "Got mutationAnnotatedSortedWindowedSum ", Sys.time())

  # x vector of quantiles representing the number of white balls drawn without replacement from an urn which contains both black and white balls.  m the number of white balls in the urn.  n the number of black balls in the urn.  k
  # the number of balls drawn from the urn.  p probability, it must be between 0 and 1.  dhyper(x, m, n, k, log = FALSE)

  # slidingWindowSizeVector <- rep(slidingWindowSize, length(mutationAnnotatedSortedWindowedSum)) if (length(mutationAnnotatedSortedWindowedSum) != length(sum(mutationAnnotatedSortedWindowedSum)) | length(slidingWindowSizeVector) !=
  # length( (length(mutationAnnotatedSortedWindowed) - sum(mutationAnnotatedSortedWindowedSum)))) # browser()

  # white balls == mutated

  # browser()
  x <- zoo::coredata(mutationAnnotatedSortedWindowedSum) # white balls draw / mutated probes every drawn
  m <- sum(mutationAnnotatedSortedLocal$MUTATIONS) # white balls total / mutated probes total
  n <- (length(mutationAnnotatedSortedLocal$MUTATIONS) - sum(mutationAnnotatedSortedLocal$MUTATIONS)) # black balls total /non mutated probes total
  k <- slidingWindowSize # number of balls drawn / number of probes drawn

  lesionpValue <- stats::dhyper(x, m, n, k)

  message(sampleName, " ", "Got lesionpValue ", Sys.time())

  #--- remove missed pValue from hypergeometric function
  lesionpValue[is.nan(lesionpValue)] <- 1

  message(sampleName, " ", "Replaces NaN  pValue in lesionpValue ", Sys.time())

  ## correction by Bonferroni
  lesionWeighted <- ((lesionpValue ) < (bonferroniThreshold/length(lesionpValue)))

  if (dim(probeFeatures)[1] != length(lesionWeighted)) {
    browser()
  }

  lesionWeighted <- data.frame(as.data.frame(probeFeatures), LESIONS = lesionWeighted)

  lesionWeighted <- sortByCHRandSTART(lesionWeighted)
  lesionWeighted <- subset(lesionWeighted, LESIONS == 1)[, c("CHR", "START", "END")]

  if (dim(lesionWeighted)[1] > dim(mutationAnnotatedSortedLocal)[1]) {
    browser()
  }

  return(lesionWeighted)
}
