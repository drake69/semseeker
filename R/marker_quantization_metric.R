#' @importFrom doRNG %dorng%
#' @importFrom doFuture %dofuture%
marker_quantization_metric <- function()
{

  # result_folder, maxResources = 90, parallel_strategy  = "multisession", ...
  # browser()
  # ssEnv <- init_env( result_folder =  result_folder, maxResources =  maxResources, parallel_strategy  =  parallel_strategy, start_fresh = FALSE, ...)
  ssEnv <- get_session_info()

  keys <- ssEnv$keys_areas_subareas_markers_figures
  keys <- keys[keys$AREA == "PROBE",]
  result_folderPivot <- dir_check_and_create(ssEnv$result_folderData,"Pivots")

  # result_folderPivot <- "~/Documents/Dati_Lavoro/osteoporosis/results/GSE99624/quantized/Pivots/"
  # key <- data.frame("MARKER" = "MUTATIONS", "FIGURE"="BOTH", "AREA"="PROBE","SUBAREA"="WHOLE")

  keys <- keys[keys$MARKER != "DELTAS",]
  keys <- keys[keys$MARKER != "DELTAR",]
  keys <- keys[keys$MARKER != "LESIONS",]
  keys <- keys[keys$MARKER != "SIGNAL",]

  keys <- keys[complete.cases(keys),]
  nkeys <- nrow(keys)

  annotate_bed()
  create_excel_pivot()

  if(ssEnv$showprogress)
    progress_bar <- progressr::progressor(along = 1:nkeys)
  else
    progress_bar <- ""

  # browser()
  to_export <- c("keys","result_folderPivot","ssEnv","file_path_build","progress_bar", "progression_index","progression","progressor_uuid","owner_session_uuid","trace")
  result_temp <- data.frame()
  # result_temp <- foreach::foreach(k = 1:nkeys, .combine =  plyr::rbind.fill, .export = to_export) %dorng%
  for (k in 1:nkeys)
  {

    # k <- 1
    key <- keys [k,]
    message(key)

    res_temp <- data.frame("MARKER" = key$MARKER, "FIGURE" = key$FIGURE, "AREA" = key$AREA, "SUBAREA" = key$SUBAREA)

    if(is.na(key$MARKER))
    {
      browser()
      next
    }

    if(key$MARKER=="DELTARP" | key$MARKER=="DELTARQ")
      original_marker <- "DELTAR"
    else
      original_marker <- "DELTAS"

    # browser()
    res_temp$ORIGINAL_MARKER <- original_marker

    pivot_subfolder <- dir_check_and_create(result_folderPivot, original_marker)
    fname <- file_path_build( pivot_subfolder ,c(original_marker, key$FIGURE, key$AREA,key$SUBAREA),"csv", add_gz=TRUE)
    if(!file.exists(fname))
    {
      log_event("DEBUG:", format(Sys.time(), "%a %b %d %X %Y")," File not found: ", fname)
      next
    }
    original <- utils::read.table(gzfile(fname), sep  =  ";", header  =  T)
    original <- original[-1, -1]
    original <- as.numeric(as.vector(as.matrix(original)))


    pivot_subfolder <- dir_check_and_create(result_folderPivot, key$MARKER)
    fname <- file_path_build( pivot_subfolder ,c(key$MARKER, key$FIGURE, key$AREA,key$SUBAREA),"csv", add_gz=TRUE)
    if(!file.exists(fname))
    {
      log_event("DEBUG:", format(Sys.time(), "%a %b %d %X %Y")," File not found: ", fname)
      next
    }
    quantized <- utils::read.table(gzfile(fname), sep  =  ";", header  =  T)
    quantized <- quantized[-1, -1]
    quantized <- as.numeric(as.vector(as.matrix(quantized)))

    mdl_perf <- model_performance(original, quantized,c(),c())
    res_temp <- cbind(res_temp, mdl_perf)

    # SSIM con pracma
    ssim_value <- ssim(original, quantized)
    res_temp$Structural_Similarity_Index <- ssim_value

    # VI
    variation_of_information <- function(original, quantized) {

      entropy <- function(p) {
        p <- p[p > 0]  # Rimuovi zeri per evitare log(0)
        -sum(p * log2(p))
      }

      joint_table <- table(original, quantized)
      joint_prob <- joint_table / sum(joint_table)
      original_prob <- margin.table(joint_prob, 1)
      quantized_prob <- margin.table(joint_prob, 2)

      H_original <- entropy(original_prob)
      H_quantized <- entropy(quantized_prob)
      H_joint <- entropy(joint_prob)

      VI <- H_original + H_quantized - 2 * H_joint
      return(VI)
    }

    vi_value <- variation_of_information(original, quantized)
    res_temp$variation_of_information <- vi_value

    # calculate JSD
    # Combine the unique elements from both samples to create a common event space
    common_events <- unique(c(original, quantized))

    # Create adjusted frequency tables for both samples
    frequency_table1_adjusted <- tabulate(match(original, common_events), nbins = length(common_events))
    frequency_table2_adjusted <- tabulate(match(quantized, common_events), nbins = length(common_events))

    # Convert adjusted frequencies to probabilities
    probability_distribution1_adjusted <- frequency_table1_adjusted / sum(frequency_table1_adjusted)
    probability_distribution2_adjusted <- frequency_table2_adjusted / sum(frequency_table2_adjusted)

    # Calculate the Jensen-Shannon distance
    jsd <- suppressMessages(suppressWarnings(philentropy::JSD(rbind(probability_distribution1_adjusted, probability_distribution2_adjusted))))
    res_temp$JSD <- jsd

    if(ssEnv$showprogress)
      progress_bar(sprintf("Doing comparison."))

    # res_temp
    result_temp <- plyr::rbind.fill(result_temp, res_temp)

    # Errore Assoluto Medio (MAE)
    # •	Willmott, C. J., & Matsuura, K. (2005). Advantages of the mean absolute error (MAE) over the root mean square error (RMSE) in assessing average model performance. Climate Research, 30(1), 79-82.
    # 2.	Indice di Similarità Strutturale (SSIM)
    # •	Wang, Z., Bovik, A. C., Sheikh, H. R., & Simoncelli, E. P. (2004). Image quality assessment: From error visibility to structural similarity. IEEE Transactions on Image Processing, 13(4), 600-612.
    # •	Avanaki, M. R. N. (2009). Exact calculation of the structural similarity index for image quality assessment. Journal of Mathematical Imaging and Vision, 34(2), 121-131.
    # 3.	Variazione dell’Informazione (VI)
    # •	Meilă, M. (2007). Comparing clusterings—an information based distance. Journal of Multivariate Analysis, 98(5), 873-895.
    # •	Vinh, N. X., Epps, J., & Bailey, J. (2010). Information theoretic measures for clusterings comparison: Variants, properties, normalization and correction for chance. Journal of Machine Learning Research, 11(10), 2837-2854.


  }
  browser
  colnames(result_temp) <- toupper(colnames(result_temp))
  dataFolder <- dir_check_and_create(ssEnv$result_folderData,c("Distributions"))
  filename  =  file_path_build(dataFolder,c("DISTRIBUTION", "ANALYSIS"),"csv")
  utils::write.csv(result_temp, file = filename, row.names = FALSE)
}


# Funzione per calcolare l'indice di similarità strutturale (SSIM) tra due vettori
ssim <- function(x, y, K1 = 0.01, K2 = 0.03, L = 1) {
  # Assicurati che i vettori siano della stessa lunghezza
  if (length(x) != length(y)) {
    stop("I vettori devono avere la stessa lunghezza.")
  }

  # Calcola le medie
  mu_x <- mean(x)
  mu_y <- mean(y)

  # Calcola le varianze
  sigma_x2 <- var(x)
  sigma_y2 <- var(y)

  # Calcola la covarianza
  sigma_xy <- cov(x, y)

  # Costanti per stabilizzare la divisione con denominatori piccoli
  C1 <- (K1 * L)^2
  C2 <- (K2 * L)^2

  # Calcolo dei tre termini del SSIM
  luminance <- (2 * mu_x * mu_y + C1) / (mu_x^2 + mu_y^2 + C1)
  contrast <- (2 * sqrt(sigma_x2) * sqrt(sigma_y2) + C2) / (sigma_x2 + sigma_y2 + C2)
  structure <- (sigma_xy + C2 / 2) / (sqrt(sigma_x2) * sqrt(sigma_y2) + C2 / 2)

  # Calcola SSIM
  ssim_value <- luminance * contrast * structure
  return(ssim_value)
}
