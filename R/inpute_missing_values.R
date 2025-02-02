inpute_missing_values <- function(signal_data){

  ssEnv <- get_session_info()
  # count number of na
  n_na <- sum(is.na(signal_data))

  n_item = nrow(signal_data)*ncol(signal_data)/100
  log_event("INFO:", format(Sys.time(), "%a %b %d %X %Y") ," Imputing missing values using ", ssEnv$inpute , " method. Number of missing values: ", n_na, " corresponding to: ", round(n_na/n_item, 2), " % of the data.")
  if (ssEnv$inpute=="median")
  {
    # Assuming signal_data is a matrix or data frame
    chunk_size <- 100000  # Define a chunk size
    for (i in seq(1, nrow(signal_data), by = chunk_size)) {
      # Define the chunk
      chunk_indices <- i:min(i + chunk_size - 1, nrow(signal_data))

      # Process the chunk
      row_medians <- apply(signal_data[chunk_indices, ], 1, stats::median, na.rm = TRUE)
      signal_data[chunk_indices, ][is.na(signal_data[chunk_indices, ])] <-
        row_medians[row(signal_data[chunk_indices, ])[is.na(signal_data[chunk_indices, ])]]
    }
  }
  else if (ssEnv$inpute=="mean")
  {
    chunk_size <- 100000  # Define a chunk size
    for (i in seq(1, nrow(signal_data), by = chunk_size)) {
      # Define the chunk
      chunk_indices <- i:min(i + chunk_size - 1, nrow(signal_data))
      # Process the chunk
      row_medians <- apply(signal_data[chunk_indices, ], 1, stats::mean, na.rm = TRUE)
      signal_data[chunk_indices, ][is.na(signal_data[chunk_indices, ])] <-
        row_medians[row(signal_data[chunk_indices, ])[is.na(signal_data[chunk_indices, ])]]
    }
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
  n_na <- sum(is.na(signal_data))
  log_event("INFO:", format(Sys.time(), "%a %b %d %X %Y") ," Imputation completed. Number of missing values: ", n_na)
  return(signal_data)
}
