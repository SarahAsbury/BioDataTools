# Import libraries --------------------------------------------------------
library(tidyverse)
library(pls)



# Subfunctions ------------------------------------------------------------
# === Sub-functions ===
#extract results from pls objects (general function)
extract.comps <- function(pls.extract, #dataframe from plsr object; columns expected in following format: responsevariable.comp.x
                          ncomp, #number of components/latent variables (plsr object)
                          response.vars, #character vector of response variable column names
                          min = TRUE, #extract minimum value
                          intercept = TRUE, #logical; is responsevariable.Intercept. a column in pls.extract?
                          transpose = TRUE
                          ){
  #inspect data
  print("Wrangling dataframe:")
  print(pls.extract[1:5,1:5])
  omitted.row <- nrow(pls.extract) - 5
  omitted.col <- ncol(pls.extract) - 5
  if(omitted.col > 0){
    print(paste("Print omitted", omitted.col, "columns."))
  }
  if(omitted.row > 0){
    print(paste("Print omitted", omitted.row, "rows"))
  }
  print("plsr extract dataframe column names (sample)")
  print(colnames(pls.extract)[1:(length(response.vars)*ncomp*0.25)])

  #initiate loop lists
  comp.performance <- list()
  comp.performance.dfs <- list()
  comp.performance.min <- list()
  for(i in 1: length(response.vars)){
    #dynamic variable
    response.temp <- response.vars[i]

    #all possible column names corresponding to response variable i
    if(intercept == TRUE){
      all.response.colnames <- c(paste(response.temp, ".Intercept.", sep = "."), paste(response.temp, 2:ncomp, "comps", sep = "."))
    }

    if(intercept == FALSE){
      all.response.colnames <- paste(response.temp, 2:ncomp, "comps", sep = ".")
      }

    # === df extract ===
    #extract response variable (e.g taxa)
    df.temp <- pls.extract %>%
      select(all_of(all.response.colnames)) %>%
      rename_with(~ gsub(pattern = paste0(response.temp, "."), replacement = "", .x))

    #fix .Intercept.
    if(intercept == TRUE){
      df.temp <- df.temp %>% rename(Intercept = .Intercept.)
    }

    #transpose
    if(transpose == TRUE){
      df.temp <- df.temp %>%
        t() %>% data.frame %>%
        rownames_to_column(var = "comps")
    }

    #assign df to list
    comp.performance.dfs[[i]] <- df.temp
    names(comp.performance.dfs)[[i]] <- response.temp

    # === min extract ===
    if(min == TRUE){
      #extract min
      minCV.comp <- df.temp %>% slice_min(adjCV, n = 1)

      #assign min to list
      comp.performance.min[[i]] <- minCV.comp
      names(comp.performance.min)[[i]] <- response.temp


    }
    # === combine lists ===
    if(min == TRUE){
      comp.performance <- list(comp.performance.dfs, comp.performance.min)
      names(comp.performance) <- c("dfs", "minCV")
    }
    if(min == FALSE){
      comp.performance <- comp.performance.dfs
    }

  }
  return(comp.performance)
}



#extract rmse from pls
pls.model.rmse <- function(pls, #pls object
                           response.vars, #character vector of prs response variables
                           intercept = TRUE, #logical; is response variable.Intercept. a column in pls.extract?
                           ncomp  #number of components/latent variables (plsr object)
)
{
  rmsep <- RMSEP(pls)
  rmsep.extract <- rmsep$val %>% data.frame
  rmse.performance <- extract.comps(pls.extract = rmsep.extract,
                                    response.vars = response.vars,
                                    min = TRUE,
                                    intercept = intercept,
                                    ncomp = ncomp)
  return(rmse.performance)
}



#plot rmse
pls.plot.rmse <- function(rmse.performance #pls.model.rmse output object
)
{
  #wrangle plot df
  plot.df <- bind_rows(rmse.performance$dfs, .id = "column_label") %>%
    mutate(comps = gsub(".comps", "", comps)) %>%
    mutate(comps = as.numeric(gsub("Intercept", 0,  comps)))

  #plot
  p <- ggline(plot.df, x = "comps", y = "adjCV", color = "column_label") +
    facet_wrap(. ~ column_label, scales = "free_y") +
    theme(legend.position = "none",
          strip.text.x = element_text(size = 7))

  return(p)
}



#extract model accuracy from pls
pls.model.accuracy <- function(rmse.performance, #pls.model.rmse output object
                               response.mat, #response matrix input to plsr function
                               pls #pls object
){
  #best # components for each taxa
  #remove response variables (e.g taxa) where intercept (0) is best performing number of components
  best.comp <- bind_rows(rmse.performance$minCV, .id = "column_label") %>% filter(comps != "Intercept") %>%
    mutate(column.select = paste(column_label, comps, sep = "."))

  #actual taxa PRS
  actual <- taxa.mat %>% data.frame %>% select(all_of(best.comp$column_label))

  #predicted taxa PRS
  predictions <- (pls$validation$pred) %>% data.frame %>%
    select(all_of(best.comp$column.select))
  accuracy.df <- cbind(actual, predictions)

  accuracy <- list(accuracy.df, best.comp)
  names(accuracy) <- c("accuracy.df", "viable.response.comp")
  return(accuracy)
}



#plot model accuracy
pls.plot.accuracy <- function(accuracy #pls.model.accuracy output object
){
  viable.response <- (accuracy$viable.response.comp$column_label)

  if(length(viable.response) > 0){
  plot.list <- list()
  for(i in 1:length(viable.response)){
    temp.response <- viable.response[i]
    accuracy.temp <- accuracy$accuracy.df %>% select(starts_with(temp.response)) %>%
      rename_with( ~ gsub(pattern = "\\.[0-9]+\\.comps", replacement = ".pred", x = .))
    print(head(accuracy.temp))
    p <- ggscatter(data = accuracy.temp,
                   x = paste(temp.response), y = paste(paste0(temp.response, ".pred")),
                   color = "black",
                   alpha = 0.5,
                   add = "reg.line",
                   conf.int = TRUE, # Add confidence interval
                   cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                   cor.coeff.args = list(method = "pearson"),
    )
    plot.list[[i]] <- p
  }
  plots <- cowplot::plot_grid(plotlist = plot.list, ncol = 2, scale = 0.95)
  }

  else{
    plots <- NULL
    warning("No response variable met criteria for viable response (i.e. Intercept had best RMSE performance). No accuracy plots returned.")
  }
  return(plots)

}




# === plot co-efficients ===
#Single co-efficient plot
pls.coef.plot <- function(pls, #object
                          response.colname #string
)
{
  if(length(response.colname) != 1){
    stop("There should be a single response variable column name provided.")
  }

  pls.coefficients <- pls$coefficients %>% data.frame() %>%
    mutate(response = (.data[[response.colname]]) %>% as.numeric) %>%
    select(response) %>%
    rownames_to_column(var = "variable")

  #plot
  p <- pls.coefficients  %>%
    ggdotchart(y = "response", x = "variable",
               ggtheme = theme_pubr(),
               add = "segment",
               color = "response") +
    scale_colour_gradientn(colours = hcl.colors(n = 10, palette = "Blue-Red")) +
    coord_flip() +
    theme(text = element_text(size = 7)) +
    ylab(response.colname)

  return(p)
}

#Looped; plot all viable response variables
pls.coefficients <- function(accuracy, #pls.model.accuracy output object; used to determine which response variables are viable
                             pls #pls object; used to extract co-efficients
                             )
  {

  #list viable response variable
  viable.response <- accuracy$viable.response.comp$column.select

  #initialize loop variables
  plot.list <- list()

  #loop
  if(length(viable.response) > 0){
    for(i in 1:length(viable.response)){
      response.colname <- viable.response[i]
      plot.list[[i]] <- pls.coef.plot(pls = pls,
                                      response.colname = response.colname)
    }
    plots <- cowplot::plot_grid(plotlist = plot.list, ncol = 2, scale = 0.95)
  }

  else{
    plots <- NULL
    warning("No response variable met criteria for viable response (i.e. Intercept had best RMSE performance). No co-efficient plots returned.")
  }

  return(plots)
}



# Pipeline Functions ------------------------------------------------------
# === Pipeline function ===
pls.pipeline <- function(response.mat, #matrix of predictor variables
                         pred.mat, #matrix of response variables,
                         ncomp = 10, #number of latent variables to try
                         validation = "LOO", #validation methods; default is Leave One Out)
                         response.vars, #character vector of response variables
                         z_scale = TRUE, #z-score scaling, default is to scale
                         intercept = TRUE #logical; is responsevariable.Intercept. a column in pls.extract?
)
{
  # === Run PLS model ===
  print("1. Run PLS Model")
  pls <- plsr(response.mat ~ pred.mat, ncomp = ncomp, validation = validation, scale = z_scale, center = z_scale)

  # === RMSE: Extract and plot ===
  print("2. Extract model performance (RMSE")
  rmse.performance <- pls.model.rmse(pls = pls, response.vars = response.vars, intercept = intercept, ncomp = ncomp)
  print("3. Plot and store model performance figure")
  rmse.plot <- pls.plot.rmse(rmse.performance = rmse.performance)

  # === Model accuracy: Extract and plot ===
  print("4. Extract model accuracy")
  accuracy <- pls.model.accuracy(rmse.performance = rmse.performance,
                                 response.mat = response.mat,
                                 pls = pls)
  print("5. Plot and store model accuracy figure")
  accuracy.plot <- pls.plot.accuracy(accuracy = accuracy)


  # === Co-efficient Plots ===
  print("6. Plot predictor variable co-efficients")
  coefficient.plot <- pls.coefficients(accuracy = accuracy, pls = pls)

  # === return results ===
  print('7. Store and return PLS results')
  pls.output <- list(pls, rmse.performance, rmse.plot, accuracy, accuracy.plot, coefficient.plot)
  names(pls.output) <- c("pls", "rmse", "rmse.plot", "accuracy", "accuracy.plot", "coefficient.plot")
  return(pls.output)


}



