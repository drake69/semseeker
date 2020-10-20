#' Title
#'
#' @param mutationAnnotatedSorted
#' @param slidingWindowSize
#' @param sampleName
#' @param probePositions
#' @param bonferroniThreshold
#'
#' @return
#' @export
#'
#' @examples
getLesions <- function(mutationAnnotatedSorted, slidingWindowSize, sampleName, probePositions, bonferroniThreshold = 0.05) {

  colnames_required <- c("CHR", "START", "END")
  check_colnames(probePositions, colnames_required, error = TRUE)
  check_colnames(mutationAnnotatedSorted, "MUTATIONS", error = TRUE)

  mutationAnnotatedSortedLocal <- sortByCHRandSTART(mutationAnnotatedSorted)
  probePositions <- sortByCHRandSTART(probePositions)

  if (!test_match_order(row.names(mutationAnnotatedSortedLocal), probePositions$Probe))
    stop("Wrong order matching Probes and Mutation!", Sys.time())

  mutationAnnotatedSortedLocal$CHR <- droplevels(mutationAnnotatedSortedLocal$CHR)

  # if (length(unique(attributes(mutationAnnotatedSortedLocal$CHR)$levels)) != 1) browser()

  # print(unique(attributes(mutationAnnotatedSortedLocal$CHR)$levels))

  if (length(mutationAnnotatedSortedLocal$CHR) <= slidingWindowSize)
    browser()

  missedWindowLength <- ((slidingWindowSize - 1)/2)

  # browser()
  mutationAnnotatedSortedWindowed <- zoo::zoo(x = mutationAnnotatedSortedLocal[, "MUTATIONS"])
  message(sampleName, " ", "Got mutationAnnotatedSortedWindowed ", Sys.time())

  sumIntoCentreWindow <- function(x) {
    zoo::rollapply(x, width = slidingWindowSize, by = 1, FUN = sum, align = "center")
  }

  ## check burden into the window (sum)
  mutationAnnotatedSortedWindowedSum <- sumIntoCentreWindow(mutationAnnotatedSortedWindowed)

  message(sampleName, " ", "Got mutationAnnotatedSortedWindowedSum ", Sys.time())

  # x vector of quantiles representing the number of white balls drawn without replacement from an urn which contains both black and white balls.  m the number of white balls in the urn.  n the number of black balls in the urn.  k
  # the number of balls drawn from the urn.  p probability, it must be between 0 and 1.  dhyper(x, m, n, k, log = FALSE)

  # slidingWindowSizeVector <- rep(slidingWindowSize, length(mutationAnnotatedSortedWindowedSum)) if (length(mutationAnnotatedSortedWindowedSum) != length(sum(mutationAnnotatedSortedWindowedSum)) | length(slidingWindowSizeVector) !=
  # length( (length(mutationAnnotatedSortedWindowed) - sum(mutationAnnotatedSortedWindowedSum)))) browser()

  # white balls == mutated

  x <- mutationAnnotatedSortedWindowedSum  # white balls drawm / mutated probes every drawn
  m <- sum(mutationAnnotatedSortedWindowedSum)  # white balls total / mutated probes total
  n <- (length(mutationAnnotatedSortedWindowed) - sum(mutationAnnotatedSortedWindowed))  # black balls total /non mutated probes total
  k <- slidingWindowSize  # number of balls drawn / number of probes drawn

  lesionpValue <- dhyper(x, m, n, k)

  # lesionpValue <- phyper(mutationAnnotatedSortedWindowedSum, sum(mutationAnnotatedSortedWindowedSum), (length(mutationAnnotatedSortedWindowed) - sum(mutationAnnotatedSortedWindowedSum)), slidingWindowSize) lesionpValue <- round(
  # lesionpValue , 10)

  lesionpValue <- coredata(lesionpValue)

  startingValue <- lesionpValue[1:missedWindowLength]
  startingValue <- 1

  endingValue <- lesionpValue[(length(lesionpValue) - missedWindowLength):length(lesionpValue)]
  endingValue <- 1

  message(sampleName, " ", "Got lesionpValue ", Sys.time())
  #--- remove missed pValue from hypergeometric function
  lesionpValue[is.nan(lesionpValue)] <- 1
  message(sampleName, " ", "Replaces NaN  pValue in lesionpValue ", Sys.time())

  ## correction by Bonferroni
  lesionWeighted <- (lesionpValue < (bonferroniThreshold/length(lesionpValue)))

  missedValue <- rep(FALSE, (slidingWindowSize - 1)/2)
  row.names(missedValue) <- row.names(startingValue)
  lesionWeighted <- append(missedValue, lesionWeighted)
  row.names(missedValue) <- row.names(endingValue)
  lesionWeighted <- append(lesionWeighted, missedValue)

  if (dim(probePositions)[1] != length(lesionWeighted))
    browser()

  lesionWeighted <- data.frame(as.data.frame(probePositions), LESIONS = lesionWeighted)

  lesionWeighted <- sortByCHRandSTART(lesionWeighted)
  lesionWeighted <- subset(lesionWeighted, LESIONS == 1)[, c("CHR", "START", "END")]

  if (dim(lesionWeighted)[1] > dim(mutationAnnotatedSortedLocal)[1])
    browser()

  if (dim(lesionWeighted)[1] > dim(subset(mutationAnnotatedSortedLocal, MUTATIONS == 1))[1]) {

    # ext <- function(i) { # x vector of quantiles representing the number of white balls drawn without replacement from an urn which contains both black and white balls.  # m the number of white balls in the urn.  # n the number of
    # black balls in the urn.  # k the number of balls drawn from the urn.  d <- mutationAnnotatedSortedWindowedSum >= i x <- as.numeric(d == 1) # white balls drawm / mutated probes every drawn m <- sum(d == 1) # white balls total /
    # mutated probes total n <- (length(d) - sum(d == 1)) # black balls total /non mutated probes total (all balls less the white balls) k <- 1 # number of balls drawn / number of probes drawn lesionpValue_1 <- dhyper(x , m , n , k)
    # lesionpValue_1[is.nan(lesionpValue_1)] <- 1 lesionWeighted_1 <- (lesionpValue_1 < (0.05/length(lesionpValue_1))) print(sum(lesionWeighted_1)) } browser() ext(1) ext(2) ext(3) ext(4) ext(5) ext(6) ext(7) ext(8) ext(9)
  }

  return(lesionWeighted)
}
