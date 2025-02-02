# sample_sheet <- utils::read.csv2("~/Desktop/What_s_going_On/papers_writing/Dioxin/orange/analysis/epimutation/Data/sample_sheet_result.csv")
# mvalue = readRDS(file = "~/Desktop/What_s_going_On/papers_writing/Dioxin/orange/analysis/champ_junone/normalizedData.rds")
# mvalue <- as.data.frame(mvalue)
# sample_sheet <- sample_sheet[,c("Sample_ID","Sample_Group","age.ch1","bmi.ch1","smoking.lvl")]
# colnames(sample_sheet) <- c("Sample_ID","Sample_Group","age.ch1","bmi.ch1","smoking.lvl")
# signal_values <-as.data.frame(mvalue)
# covariates <- c("age.ch1","bmi.ch1","smoking.lvl")
# library(doRNG)
#' @importFrom doRNG %dorng%
covariates_removal <- function(signal_values, covariates, sample_sheet, result_folder, ...)
{
  init_env( result_folder= result_folder, maxResources= maxResources, ...)
  ssEnv <- get_session_info()
  utils::write.csv2(sample_sheet, file = file.path(ssEnv$result_folderData, "sample_sheet.csv"))
  # sample_sheet <- subset(sample_sheet, Sample_Group != "Reference")
  signal_values <- signal_values[,sample_sheet$Sample_ID]
  get_th <- function(signal_row)
  {
    #
    # signal_row <- signal_values[5,]
    b_values <- as.vector(t(signal_row))
    temp_df <- cbind(b_values, sample_sheet)
    model_formula <- paste("b_values ~ ", paste(covariates, collapse = " + "))
    model <- stats::glm( formula= model_formula, data = temp_df)
    # get the coefficients
    coeffs <- summary(model)$coefficients
    coeffs <- coeffs[,1]
    p_value <- summary(model)$coefficients[,4]
    mask <- p_value < as.numeric(ssEnv$alpha)
    n_values <- length(covariates) + 1
    #
    mask <- mask[2:n_values]
    intercept <- coeffs[1] * mask[1]
    # log_event(sum(mask))
    temp_result <- signal_row
    if (sum(mask) != 0)
    {
      # calculate product of columns and mask column by column
      temp_result <- replace(temp_df[, covariates], !mask, 0)
      temp_result <- as.data.frame(temp_result)
      temp_result <- apply(temp_result * coeffs[2:n_values],1,sum)
      temp_result <- b_values - temp_result - intercept
      temp_result <- as.data.frame(t(temp_result))
    }
    temp_result <- as.data.frame(temp_result)
    colnames(temp_result) <- colnames(signal_row)
    rownames(temp_result) <- rownames(signal_row)
    temp_result
  }

  if(ssEnv$showprogress)
    progress_bar <- progressr::progressor(along =1:nrow(signal_values))
  else
    progress_bar <- ""

  to_export <- c("get_th","signal_values","covariates","sample_sheet")
  th <- foreach::foreach(i = 1:nrow(signal_values), .combine = rbind, .export = to_export) %dorng%
    # for (i in 1:nrow(signal_values))
    {
      signal_row <- signal_values[i,]
      # every 10% of the probes, log_event the progress
      ten_perc <- round(nrows(signal_values)/10)
      if(ssEnv$showprogress)
      {
        progress_bar(sprintf("Done probes: %s", i))
      }
      tmp <- get_th(signal_row)
      # if (i == 1)
      #   th <- tmp
      # else
      #   th <- rbind(th, tmp)
      tmp
    }
  saveRDS(th, file = file.path(ssEnv$result_folderData, "signal_filtered.rds"))
}
# covariates_filtering(signal_values, covariates, sample_sheet)

