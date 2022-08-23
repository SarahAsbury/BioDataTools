#General functions used for cleaning datasets


# Pipe Status Functions ----------------------------------------------------------
#pipe_underpass
#https://stackoverflow.com/questions/46123285/printing-intermediate-results-without-breaking-pipeline-in-tidyverse
pipe_underpass <- function(df, fun){
  fun(df)
  return(df)
}

#pipe_message
pipe_message <- function(df, statement){
  message(statement)
  return(df)
}

#pipe_data_status
pipe_data_status <- function(df){
  print(head(df))
  return(df)
}
