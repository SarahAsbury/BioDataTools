#General functions used for cleaning datasets


#string_as_variable
#Uses character string assigned to an R object, and treats it as a variable
#Useful for creating dynamic variables inside functions
  #Especially ggplot2

strings_as_variable <- function(x){
  return(eval(parse(text = x)))
}
