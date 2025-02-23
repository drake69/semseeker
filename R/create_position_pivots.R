#' create pivot from samples' single bed for the PROBE (WHOLE)
#'
#' @param population population data frame to use for creating pivot
#' @param keys markers and figure to create pivot of
#'
create_position_pivots <- function(population, keys)
{
  # create pivot for basic markers and figure
  ssEnv <- get_session_info()
  for (k in 1:nrow(keys))
  {
    key <- keys[k,]
    for ( s in 1:nrow(population))
    {
      sample <- population[s,]
      create_position_pivot_from_single_bed(sample$Sample_Group,sample$Sample_ID, as.character(key$MARKER), as.character(key$FIGURE), "POSITION","")
    }
  }
}
