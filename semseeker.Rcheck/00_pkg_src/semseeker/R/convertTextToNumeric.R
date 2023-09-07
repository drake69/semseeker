convertTextToNumeric <- function(text) {

  # Get the system's locale settings
  locale_settings <- Sys.localeconv()

  # Extract the decimal separator from the locale settings
  decimal_separator <- paste("[",locale_settings["decimal_point"],"]", sep="")

  # Example usage:
  # print(paste("Decimal separator:", decimal_separator))

  if(length(gregexpr("\\.", text)[[1]])>1 | length(gregexpr("\\,", text)[[1]])>1)
  {
    # Check if the text contains a decimal separator
    if (grepl("[.,]", text)) {

      # Try to convert with "." as the decimal separator
      # and  "," as thousands separator
      converted_value <- suppressWarnings(gsub("\\.", "", text))

      # If conversion successful, return the numeric value
      if (!is.na(converted_value))
        return(converted_value)

      # Try to convert with "," as the decimal separator
      converted_value <- suppressWarnings(gsub(",", "", text))

      # If conversion successful, return the numeric value
      if (!is.na(converted_value))
        return(converted_value)
    }
  }
  else
    {
      # Try to convert with "." as the decimal separator
      # and  "," as thousands separator
      converted_value <- suppressWarnings(as.numeric(gsub(",", ".", text)))

      # If conversion successful, return the numeric value
      if (!is.na(converted_value))
        return(converted_value)

      # Try to convert with "," as the decimal separator
      converted_value <- suppressWarnings(as.numeric(gsub("\\.", ",", text)))

      # If conversion successful, return the numeric value
      if (!is.na(converted_value))
        return(converted_value)
    }

  # If no decimal separator found or conversion unsuccessful, return NA
  return(NA)
}



text_value <-  "1,200,000.34"
 semseeker:::convertTextToNumeric(text_value)==1200000.34

text_value <-  "1.200.000,34"
semseeker:::convertTextToNumeric(text_value)==1200000.34
