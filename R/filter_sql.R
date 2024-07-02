filter_sql <- function(sql_conditions, data_frame)
{

  # browser()

  # check if sql_conditions is am empty string
  if (length(sql_conditions) == 0)
    return(data_frame)

  for ( s in 1:length(sql_conditions))
  {
    sql_condition <- sql_conditions[s]
    if(!is.null(sql_condition))
      if (sql_condition != "")
      {
        # sql_condition <- toupper(sql_condition)
        sql_condition <- gsub("  ", " ", sql_condition)
        if (grepl("FROM", sql_condition))
        {
          if (!grepl("FROM TABLE", sql_condition))
          {
            log_event("ERROR: ",format(Sys.time(), "%a %b %d %X %Y")," the only allowed table placeholder with sub statement is TABLE!")
            stop()
          }
          else
            sql_condition <- gsub("TABLE","data_frame", sql_condition)
        }
        sql <- paste("select * from data_frame where ", sql_condition)
        log_event("DEBUG:", format(Sys.time(), "%a %b %d %X %Y"), " Executed sql: " , sql_condition)
        columns <- toupper(colnames(data_frame))
        # drop duplicates columns
        data_frame <- data_frame[,!duplicated(columns)]
        data_frame <- sqldf::sqldf(sql)
      }
  }
  return(data_frame)
}
