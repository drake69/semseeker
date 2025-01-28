inpute_missing_values <- function(signal_data){

  ssEnv <- get_session_info()

  log_event("INFO:", format(Sys.time(), "%a %b %d %X %Y") ," Imputing missing values using ", ssEnv$inpute)
  if (ssEnv$inpute=="median")
  {
    # Assuming signal_data is a matrix or data frame
    # signal_data <- signal_data %>%
    #   dplyr::rowwise() %>%
    #   dplyr::mutate(dplyr::across(dplyr::everything(), ~ ifelse(is.na(.), median(dplyr::c_across(dplyr::everything()), na.rm = TRUE), .))) %>%
    #   dplyr::ungroup()
    chunk_size <- 10000  # Define a chunk size
    for (i in seq(1, nrow(signal_data), by = chunk_size)) {
      # Define the chunk
      chunk_indices <- i:min(i + chunk_size - 1, nrow(signal_data))

      # Process the chunk
      row_medians <- apply(signal_data[chunk_indices, ], 1, median, na.rm = TRUE)
      signal_data[chunk_indices, ][is.na(signal_data[chunk_indices, ])] <-
        row_medians[row(signal_data[chunk_indices, ])[is.na(signal_data[chunk_indices, ])]]
    }
  }
  else if (ssEnv$inpute=="mean")
  {
    signal_data <- signal_data %>%
      dplyr::rowwise() %>%
      dplyr::mutate(dplyr::across(dplyr::everything(), ~ ifelse(is.na(.), mean(dplyr::c_across(dplyr::everything()), na.rm = TRUE), .))) %>%
      dplyr::ungroup()
  }
  else if (grepl("knn", ssEnv$inpute))
  {
    if (len(strsplit(ssEnv$inpute, ";")[[1]]) != 3)
    {
      log_event("ERROR:", format(Sys.time(), "%a %b %d %X %Y") ," Invalid inpute value. Please provide the number of centers and the number of clusters separated by a semicolon.")
      stop()
    }
    centers <- strsplit(ssEnv$inpute, ";")[[1]][2]
    k <- strsplit(ssEnv$inpute, ";")[[1]][3]
    signal_matrix <- as.matrix(signal_data)
    imputed_matrix <- KMEANS.KNN::kmeans_knn(signal_matrix, centers = centers, k = k)
  }
  else if (ssEnv$inpute=="none")
  {
    signal_data <- signal_data
  }
  else
  {
    log_event("ERROR:", format(Sys.time(), "%a %b %d %X %Y") ," Invalid inpute value. Please provide one of the following: median, mean, knn.")
    stop()
  }

  gc()
  log_event("INFO:", format(Sys.time(), "%a %b %d %X %Y") ," Imputation completed.")
  return(signal_data)
}
