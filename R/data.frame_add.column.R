data.frame_add.column <- function(df,col_name, value)
{
  if(any(col_name %in% colnames(df)))
    df[,col_name] <- value
  else
  {
    tmp <- data.frame("TMP"="")
    tmp$TMP <- value
    colnames(tmp) <- col_name
    if (nrow(df)==0)
      df <- tmp
    else
      df <- cbind(df,tmp)
  }
  return(df)
}
