set_env_variable <- function(ssEnv, arguments, var_name, var_default_value)
{

  if(is.null(ssEnv[[var_name]]))
    ssEnv[[var_name]]  <- var_default_value

  if(!is.null(arguments[[var_name]]))
    ssEnv[[var_name]]  <- arguments[[var_name]]

  update_session_info(ssEnv)
  log_event("INFO: ",format(Sys.time(), "%a %b %d %X %Y"), " " , var_name ," ",  ssEnv[[var_name]])

  # remove from arguments var_name
  arguments[[var_name]] <- NULL

  return(arguments)
}
