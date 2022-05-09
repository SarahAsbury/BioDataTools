
# Load packages -----------------------------------------------------------
library(logr)
library(lubridate)


# Function: Random Forest Wrapper -----------------------------------------
#Wrapper for both regression and classification random forest
#Steps:
#1. Split test set
#2. Setup mtry parameters
#3. Run random forest
#4. Variable importance plots
#5. Top20 density loop
#6. Aggregate cm results
#7. Mean performance indicators (from best.performance)

#' @title Random Forest Wrapper
#' @description Convenient wrapper function that runs Random Forest pipeline. This pipeline include 5-fold cross-validation and hyper-parameter tuning for mtry and ntree. Number of train/test splits and train/validation ratios are customizable by user.
#' @param rf.type One of: \code{"class"} for classification trees or \code{"reg"} for regression trees.
#' @param vpred String. Name of response variable.
#' @param df Input dataframe.
#' @param dir Directory path for export. String.
#' @param rf.param Save names for folders. Input 1-word or phrase that describes: c(dataset, predictor variable)  , Default: c(dataframe.name, predictors.name)
#' @param nsets (Optional) Number of train/validation sets to generate, Default: 10
#' @param split.param (Optional) Proportion of sample assigned to training set. Range from 0 to 1. E.g 0.8 indicates 80% of sample assigned to training set for a 80:20 train:test split, Default: c(train.ratio = 0.8)
#' @param mtry (Optional) Range of mtry values to try for Random Forest hyperparameter tuning. If NA, will use mtry.guide to select optimal mtry based on tree type and number of predictor variables, Default: NA
#' @param ntree (Optional) Number of trees to try during Random Forest hyperparameter tuning, Default: (1:10) * 500
#' @param top.variables (Optional) Number of top important predictor variables to plot, Default: NA
#' @param varimp.param (Optional) Parameters to pass to variable importance plot. Default: c(selection_type = NA, metric = NA, xlab = NA)
#' @param density.param (Optional) Parameters to pass to density plot, Default: c(scale = 0.8, ncol = 2, tsize = 10, xlab = NA)
#' @param extract.names.df (Optional) Dataframe. A df must be provided if \code{density} param xlab = "extract". 1st columns is how the predictor variables should appear in the density plots (e.g IL-6 Concentration). 2nd column is how the predictor variable appears in input \code{df} (e.g IL.6.Concentration), Default: NA
#' @param experiment.note (Optional) User input human-readable note that will be sent to output log. May be used to log why/what is being run., Default: NA
#' @return Exports random forest results to sub-folders within \code{dir}
#' @details DETAILS
#' @examples
#' rf_standard(rf.type = "class", vpred = "Genotype", df = blood %>% select(-c(Sex, AnimalID)), dir = "/Users/Documents/experiment", experiment.note = "Predict mouse genotype from immune populations. No genotype excluded from dataframe. Exclude sex metadata.", rf.param = c(dataframe.name = "blood-allgenotypes", predictors.name = "immune"))
#' @rdname rf_standard
#' @export
rf_standard <- function(rf.type, #"class" or "reg",
                        vpred,
                        df,
                        dir, #parent directory to store results in
                        nsets = 10,
                        split.param = c(train.ratio = 0.80),
                        mtry = NA,
                        ntree = (1:10)*500,
                        top.variables = NA, #number of top predictor variables to plot
                        rf.param = c(dataframe.name, #name of input dataframe
                                     predictors.name), #broad 1-word description(s) of predictor variables
                        varimp.param = c(selection_type = NA, #top or random_top
                                         metric = NA,
                                         xlab = NA),
                        density.param = c(scale = 0.8, ncol = 2, tsize = 10, xlab = NA),
                        extract.names.df = NA, #provide df if xlab = "extract". 1st col is desired x label names. 2nd col is how columns appear in df input.
                        experiment.note = NA
                        )
{
  print(rf.param)
  # Setup wd  -------------------------------------------------------------------
  #Save original directory
  og.dir <- getwd() #initial directory. set back to once function done.

  #Set parent directory
  setwd(dir)

  #Get the date
  date <- lubridate::today()

  #Create and set rf sub-directory
  save.name <- paste(rf.param["dataframe.name"], date, rf.param["predictors.name"], vpred, sep = "_")
  dir.create(save.name)
  directory <- paste(dir, save.name, sep = "/")
  setwd(directory)



  # Setup log input ---------------------------------------------------------
  #Start log
  log.results <- log_open(paste0(directory, "/rf_log.txt"))

  #Date and notes
  if(!(is.na(experiment.note))){
    log_print("Experiment note:")
    log_print(experiment.note)
  }

  #Input df
  df.input <- df
  log_print("Input dataframe summary:")
  log_print(summary(df))
  log_print(df %>% count(!!as.name(vpred)))

  #Log response and predictor variables
  log_print("Response variable:")
  log_print(vpred)

  log_print("Predictor variables")
  log_print(df %>% select(-vpred) %>% colnames)



  # Set predictor variable to factor ----------------------------------------
  if(rf.type == "class"){
    df <- df %>% mutate(!!as.name(vpred) := as.factor(!!as.name(vpred)))

    #Log predictor variable set to factor
    log_print("Predictor variable set to factor")
    log_print(df %>% glimpse)
}


  # Split dataframe ---------------------------------------------------------
  log_print(
  split.na.removed <- split.ratio(df = df,
                                  nsets = nsets,
                                  train.ratio = split.param["train.ratio"],
                                  rf.type = rf.type,
                                  vpred = vpred)
  )
  #Logging: NA-removed
  log_print("Input dataset with incomplete samples removed. This is the dataset splits are perfomed on")
  log_print(split.na.removed)
  log_print(split.na.removed %>% summary)
  log_print(paste("Number of samples:", nrow(split.na.removed)))

  #Logging: split sets
  for(i in 1:nsets){
    #Set #
    log_print(paste("Set", i))

    #Get train/hold sets
    tname <- paste("train", i, sep = ".")
    train <- get(tname)
    hname <- paste("hold", i, sep = ".")
    hold <- get(hname)

    #Log
    log_print(paste("training set:", train %>% count(!!as.name(vpred))))
    log_print(paste("validation set:", hold %>% count(!!as.name(vpred))))
  }



  # Run rf ------------------------------------------------------------------
  #Set mtry parameters
  if(is.na(mtry)){
    log_print(
    mtry.guide.save <- mtry.guide(npred = df %>% select(-vpred) %>% ncol,
                                   rf.type = rf.type)
    )

    mtry <- mtry.guide.save$range
    log_print("Number of variable available at each split (mtry) automatically set to:")
    log_print(mtry)
  }

    #Log rf parameters
  log_print("Rf parameters")
  log_print(paste("Predictor variable:", vpred))
  log_print(paste("Number of train/test sets", nsets))
  log_print(paste("Hyperparameter tuning - test number of variable available at each split:"))
  log_print(mtry)
  log_print(paste("Hyperparameter tuning - test number of boostrapped aggregated trees:"))
  log_print(ntree)

  #Run rf
   #start log
  if(rf.type == "class"){
    log_print(
    rf.results <- sa_rf(vpred,
                        mtry = mtry,
                        ntree = ntree,
                        nset = nsets,
                        dataframe.name = rf.param["dataframe.name"],
                        predictors.name = rf.param["predictors.name"],
                        custom.name = TRUE,
                        date = date,
                        wd = dir
                        )
    )
  }

  if(rf.type == "reg"){
    log_print(
      rf.results <- sa_rfreg(vpred,
                             mtry = mtry,
                             ntree = ntree,
                             nset = nsets,
                             dataframe.name = rf.param["dataframe.name"],
                             predictors.name = rf.param["predictors.name"],
                             custom.name = TRUE,
                             date = date,
                             wd = dir
                             )
    )
  }



  # Re-set wd ---------------------------------------------------------------
  setwd(directory)

  # Automatically set top random variables and selection type ---------------
  log_print("Variable importance plots")
  if(is.na(top.variables)){
    top.variables <- ncol(df %>% select(-vpred))
  }

  if(top.variables > ncol(df %>% select(-vpred)))
    {
    log_print("Warning: the number of top predictor variables to plot exceeds the number of predictor variables.
              Automatically adjusting top variables to plot to be equal to the number of predictor variables.")
    top.variables <- ncol(df %>% select(-vpred))
  }

  if(is.na(varimp.param["selection_type"])){
      if(top.variables == ncol(df %>% select(-vpred))){
        varimp.param["selection_type"] <- "top"
      } else{
        varimp.param["selection_type"] <- "random_top"
      }
    }

  if(is.na(varimp.param["metric"])){
    varimp.param["metric"] <- "gini"
  }


  #Logging:
  log_print("Variable importance plots parameters:")
  log_print(paste0("Variable importance plot type:", varimp.param["selection_type"]))
  log_print(paste0("Number of top predictor variables to plot:", top.variables))
  log_print(paste0("Manually supplied x labels:", varimp.param["xlab"]))
  log_print(paste0("Metric used to rank variable importance for density plot order:", varimp.param["metric"]))
  if(!(is.na(varimp.param["xlab"]))){
    if(varimp.param["xlab"] == "extract"){
    log_print("X label name extraction dataframe:")
    log_print(extract.names.df)
     }}

  # Variable Importance Plots -----------------------------------------------
  #Automatic assignment of selection type
    #If there are more
  varimp <- rf.results$varimp
  log_print(varimp)

    varimp_plot <- varimport_plot(varimp = varimp,
                                  rf.type = rf.type,
                                  metric = "gini",
                                  selection_type = varimp.param["selection_type"],
                                  xlab = varimp.param["xlab"],
                                  extract.names.df = extract.names.df,
                                  top = top.variables)

    tiff("GiniVarimp_plot.tiff", res = 300, height = 1700, width = 1500)
    print(varimp_plot)
    dev.off()


    if(rf.type == "reg"){
    varimp_plot <- varimport_plot(varimp = varimp,
                                  rf.type = rf.type,
                                  metric = "mse",
                                  selection_type = varimp.param["selection_type"],
                                  xlab = varimp.param["xlab"],
                                  extract.names.df = extract.names.df,
                                  top = top.variables)

    tiff("MSEVarimp_plot.tiff", res = 300, height = 1700, width = 1500)
    print(varimp_plot)
    dev.off()
    }


    if(rf.type == "class"){
      varimp_plot <- varimport_plot(varimp = varimp,
                                    rf.type = rf.type,
                                    metric = "mda",
                                    selection_type = varimp.param["selection_type"],
                                    xlab = varimp.param["xlab"],
                                    extract.names.df = extract.names.df,
                                    top = top.variables)

      tiff("MDAVarimp_plot.tiff", res = 300, height = 1700, width = 1500)
      print(varimp_plot)
      dev.off()
    }




  # Density plots -------------------------------------------------------
  if(rf.type == "class"){
  log_print("Density plots")
  density <- multi.density.plot(df = split.na.removed,
                                varimp = varimp,
                                response = vpred,
                                top = top.variables,
                                rf.type = rf.type,
                                tsize = density.param["tsize"],
                                metric = varimp.param["metric"],
                                scale = density.param["scale"],
                                ncol = density.param["ncol"],
                                xlab = density.param["xlab"],
                                extract.names.df = extract.names.df
                                )

  tiff("density_plots.tiff", width = 1200, height = 3000, res = 300)
  print(density)
  dev.off()
  }



  # Aggregate confusion matrices -------------------------------------------
  if(rf.type == "class"){
    log_print("Aggregate confusion matrices")
    cm.rownames <- rf.results$rf.cm %>% rownames
    cm <- rf.results$rf.cm %>% as.data.frame %>% mutate(X = cm.rownames)

    log_print("Aggregate by count:")


    a.cm <- cm.aggregate(cm)

    log_print("Calculate proportion of aggregated confusion matrix:")
    p.cm <- cm.prop(a.cm)

    write.csv(file = "aggregated-cm.csv", x = a.cm, row.names = FALSE)
    write.csv(file = "proportion-cm.csv", x = p.cm, row.names = FALSE)
  }



  # Mean and Stdev of Performance Indicators ---------------------------------------------
  log_print("Mean and standard deviation of all train/validation set model performance.")
  model_performance <- mean_var.rfmodel(rf.type = rf.type,
                                        performance.table = rf.results$best.performance)
  log_print(model_performance)
  write.csv(x = model_performance, file = "mean-sd_modelsets.csv")

  log_close()
}
