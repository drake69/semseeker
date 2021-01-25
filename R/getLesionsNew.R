
# slidingWindowSize <- 11
# bonferroniThreshold <- 0.01
# grouping_column <- "GENE"
# mutationAnnotatedSortedLocal <- read.csv2("/home/lcorsaro/Documents/SEMSEEKER_TEST_BWS/mutations_annotated_sorted_gene.csv")

#' @importFrom dplyr %>%
getLesionsNew <- function(slidingWindowSize, bonferroniThreshold, grouping_column, mutationAnnotatedSorted)
{

  mutationAnnotatedSortedLocal <- mutationAnnotatedSorted
  summed <- aggregate(mutationAnnotatedSortedLocal$MUTATIONS, by = list(mutationAnnotatedSortedLocal[,grouping_column]), FUN = sum)
  colnames(summed) <- c(grouping_column,"MUTATIONS_COUNT")
  counted <- aggregate(mutationAnnotatedSortedLocal$MUTATIONS, by = list(mutationAnnotatedSortedLocal[,grouping_column]), FUN = length)
  colnames(counted) <- c(grouping_column,"PROBES_COUNT")
  mutationAnnotatedSortedLocal <- merge(mutationAnnotatedSortedLocal,summed, by = grouping_column)
  mutationAnnotatedSortedLocal <- merge(mutationAnnotatedSortedLocal,counted, by = grouping_column)
  rm(counted)
  rm(summed)

  #function to calculate rolled sum, returns a column vector
  roll<-function(x,lags){
    # if (length(x)<lags) {
    #   # tmp=c(rep(0,length(x)))
    #   for (i in 1:length(x))
    #   {
    #     i_inf <- max(0, i -  ((lags - 1) / 2) )
    #     i_sup<- min(length(x), (i+  ((lags - 1) / 2)))
    #     burden <- sum(x[ i_inf:i_sup ])
    #     x[i]  <- burden
    #   }
    # }
    # else {
    #
    missedWindowLength <- ((lags - 1) / 2)
    y <- c(rep(0,missedWindowLength))
    tmp <- c(y,x,y)
    tmp=zoo::rollsum( tmp, lags, align = "center", fill = 0)
    tmp <- tmp[(missedWindowLength+1):(missedWindowLength+ length(x))]
    # }
    tmp=as.numeric(tmp)
    return(tmp)
  }

  #eval(parse(text = grouping_colum))

  mutationAnnotatedSortedLocal <- mutationAnnotatedSortedLocal %>% dplyr::group_by(eval(parse(text = grouping_column))) %>% dplyr::mutate(ENRICHMENT = ave(MUTATIONS, eval(parse(text = grouping_column)), FUN = function(x) roll(x, slidingWindowSize))) %>% dplyr::ungroup ()

  mutationAnnotatedSortedLocal$ENRICHMENT[ is.na(mutationAnnotatedSortedLocal$ENRICHMENT)] <- 0

  lesionpValue <- stats::dhyper(mutationAnnotatedSortedLocal$ENRICHMENT, mutationAnnotatedSortedLocal$MUTATIONS_COUNT, mutationAnnotatedSortedLocal$PROBES_COUNT, slidingWindowSize)

  lesionpValue[is.nan(lesionpValue)] <- 1
  lesionpValue[is.na(lesionpValue)] <- 1

  tt <- data.frame(mutationAnnotatedSortedLocal,lesionpValue)

  ## correction by Bonferroni
  lesionWeighted <- ((tt$lesionpValue ) < (bonferroniThreshold/length(tt$PROBES_COUNT)))
  table(lesionWeighted)
  rm(tt)

  lesionWeighted <- data.frame(as.data.frame(mutationAnnotatedSortedLocal), LESIONS = lesionWeighted)

  lesionWeighted <- sortByCHRandSTART(lesionWeighted)
  lesionWeighted <- subset(lesionWeighted, LESIONS == 1)[, c("CHR", "START", "END")]

  if (dim(lesionWeighted)[1] > dim(mutationAnnotatedSortedLocal)[1]) {

  }

  return(lesionWeighted)

}
