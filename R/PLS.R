# Import libraries --------------------------------------------------------
library(tidyverse)
library(pls)



# Functions: PLS ----------------------------------------------------------
# === Sub-functions ===
#extract results from pls objects (general function)
extract.comps <- function(pls.extract, #dataframe from plsr object; columns expected in following format: responsevariable.comp.x
                          response.vars = taxa, #character vector of response variable column names
                          min = TRUE, #extract minimum value
                          intercept = TRUE, #logical; is responsevariable.Intercept. a column in pls.extract?
                          transpose = TRUE){
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

  #initiate loop lists
  comp.performance <- list()
  comp.performance.dfs <- list()
  comp.performance.min <- list()
  for(i in 1: length(response.vars)){
    #dynamic variable
    response.temp <- response.vars[i]

    # === df extract ===
    #extract taxa
    df.temp <- pls.extract %>%
      select(starts_with(paste(response.temp))) %>%
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
                           response.vars #character vector of prs response variables
)
{
  rmsep <- RMSEP(pls)
  rmsep.extract <- rmsep$val %>% data.frame
  rmse.performance <- extract.comps(pls.extract = rmsep.extract,
                                    response.vars = taxa,
                                    min = TRUE)
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
  #remove taxa where intercept (0) is best performing number of components
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
  plot.list <- list()
  for(i in 1:length(viable.response)){
    temp.response <- viable.response[i]
    accuracy.temp <- accuracy$accuracy.df %>% select(starts_with(temp.response)) %>%
      rename_with( ~ gsub(pattern = "[0-9]\\.comps", replacement = "pred", x = .))
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
  return(plots)
}


#plot co-efficients



# === Pipeline function ===
pls.pipeline <- function(response.mat, #matrix of predictor variables
                         pred.mat, #matrix of response variables,
                         ncomp = 10, #number of latent variables to try
                         validation = "LOO", #validation methods; default is Leave One Out)
                         response.vars, #character vector of response variables
                         z_scale = TRUE #z-score scaling, default is to scale
)
{
  # === Run PLS model ===
  pls <- plsr(response.mat ~ pred.mat, ncomp = ncomp, validation = validation, scale = z_scale, center = z_scale)

  # === RMSE: Extract and plot ===
  rmse.performance <- pls.model.rmse(pls = pls, response.vars = response.vars)
  rmse.plot <- pls.plot.rmse(rmse.performance = rmse.performance)

  # === Model accuracy: Extract and plot ===
  accuracy <- pls.model.accuracy(rmse.performance = rmse.performance,
                                 response.mat = response.mat,
                                 pls = pls)
  accuracy.plot <- pls.plot.accuracy(accuracy = accuracy)

  # === return results ===
  pls.output <- list(pls, rmse.performance, rmse.plot, accuracy, accuracy.plot)
  names(pls.output) <- c("pls", "rmse", "rmse.plot", "accuracy", "accuracy.plot")
  return(pls.output)
}



