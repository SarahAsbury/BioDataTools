#Functions to run random forest and process results



# Load Packages -----------------------------------------------------------
library(tidyverse)
library(ggpubr)
library(cowplot)
library(randomForest)
library(ViewPipeSteps) #    %P>% = print and pipe operator
library(e1071)

# Functions: supporting rf ---------------------------------------------------
#' @title Recommended mtry range
#' @description Recommended mtry range based on number of paramerets and tree type
#' @param npred Numeric. Input number of predictor variables.
#' @param rf.type One of: \code{"class} for classification trees or \code{"reg"} for regression trees.
#' @param spread.upper Upper range of mtry to use above optimal mtry value, Default: 5
#' @param spread.lower Lower range of mtry to use below optimal mtry value: 4
#' @return Returns a list: \cr
#'  \code{ideal.mtry} is optimal mtry based on tree type and number of predictor variables \cr
#'  \code{upper.mtry} is upper range of mtry for downstream hyperparameter tuning \cr
#'  \code{lower.mtry} is lowest value of mtry for downstream hyperparameter tuning \cr
#'  \code{range} is range of mtry values to try in downstream hyperparameter tuning given in the format: \code{upper.mtry:lower.mtry}
#'
#'
#' @details
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname mtry.guide
#' @export
mtry.guide <- function(npred, #Numeric. Number of predictor variables.
                       rf.type, #Character. One of: reg, class. Whether using regression (reg) or classification (class).
                       spread.upper = 5, #Total number of mtry to attempt above ideal mtry in hyperparameter tuning
                       spread.lower = 4 #Total number of mtry to attempt below ideal mtry in hyperparameter tuning
)
{
  if(rf.type == "class"){
    m <- sqrt(npred) #sqrt(p) for mtry classification guide
  }

  if(rf.type == "reg"){
    m <- npred/3 #p/3 for mtry regression guide
  }

  m.upper <- round(m + spread.upper, digits = 0)
  m.lower <- round(m - spread.lower, digits = 0)

  #Unit test: m.lower greater than 0?
  if(m.lower <= 0){
    print("Warning: lower range of mtry is less than 0. Adjusting mtry range all values are greater than 0. Consider reducing the mtry guide spread.")
    m.add <- 1 - m.lower
    m.lower <- m.lower + m.add
    m.upper <- m.upper + m.add

    print(paste0("Added: ", m.add))
  }

  range <- eval(parse(text = paste0(m.lower, ":", m.upper)))

  return(list(ideal.mtry = m, upper.mtry = m.upper, lower.mtry = m.lower, range = range))
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param rf.cm PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname rf.resplot
#' @export
rf.resplot <- function(rf.cm)
  #results plot for random forest regression
{
  accuracy.cor <- cor(rf.cm$actual, rf.cm$pred)
  cor.sign <- ifelse(accuracy.cor > 0, "pos", "neg")
  print(accuracy.cor)
  print(cor.sign)

  resplot <- rf.cm %>%
    ggscatter(x = "actual", y = "pred",
              add = "reg.line", conf.int = TRUE,
              add.params = list(color = "blue",
                                fill = "lightgray")) +
    stat_cor(method = "pearson", aes(label = ..r.label..),
             label.x = max(rf.cm$actual) * 0.8,
             label.y = ifelse(cor.sign == "pos",
                              max(rf.cm$pred) * 0.15,
                              max(rf.cm$pred) * 0.85)
    ) +
    geom_abline(intercept = 0, slope = 1, linetype = 2) +
    xlim(0, max(rf.cm$actual)) +
    ylim(0, max(rf.cm$pred))
}


# Functions: Split Methods ------------------------------------------------
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param df PARAM_DESCRIPTION
#' @param nsets PARAM_DESCRIPTION, Default: 10
#' @param train.ratio PARAM_DESCRIPTION, Default: 0.8
#' @param rf.type PARAM_DESCRIPTION
#' @param vpred PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname split.ratio
#' @export
split.ratio <- function(df,
                        nsets = 10, #number of train/tests
                        train.ratio = 0.80, #proportion; default is 80:20
                        rf.type, #One of: "class" or "reg"
                        vpred #Intended predictor variable. If rf.type = class, this is a required input.
)
{
  #Unit test: ratio
  if(train.ratio > 1){
    print("Error: Train.ratio must be a proportion")
    print("Running UNDECLARED to kill script.")
    UNDECLARED()
  }

  #Unit test: check that predictor variable is a factor
  if(rf.type == "class"){
    qc <- is.factor(eval(parse(text = paste0("df$", vpred))))
    if(qc == FALSE){
      print("Error: Predictor variable must be a factor.")
      print("Running UNDECLARED to kill script.")
      UNDECLARED()
    }
  }

  #Split
  input <- 1:nsets
  tsets <- as.character()
  hsets <- as.character()
  x2 <- df
  print(sapply(x2, function(x) sum(is.na(x2))))
  x3 <- df %>% na.omit
  print(sapply(x3, function(x) sum(is.na(x3))))


  for (i in input){
    set.seed(i)
    print(i)

    sample <- sample(1:nrow(x3), nrow(x3)*train.ratio)
    assign(paste("train", i, sep = "."), x3[sample,], envir = .GlobalEnv) #train.x (using tune wrapper downstream to separate into train vs. validation)
    assign(paste("hold", i, sep = "."), x3[-sample,], envir = .GlobalEnv) #hold.x; test/hold-out set
  }

  return(x3)

}





# Function: Regression RF ----------------------------------------------------
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param vpred PARAM_DESCRIPTION
#' @param custom.name PARAM_DESCRIPTION, Default: FALSE
#' @param date PARAM_DESCRIPTION, Default: ''
#' @param dataframe.name PARAM_DESCRIPTION, Default: ''
#' @param predictors.name PARAM_DESCRIPTION, Default: ''
#' @param mtry PARAM_DESCRIPTION, Default: 1:10
#' @param ntree PARAM_DESCRIPTION, Default: (1:10) * 500
#' @param nset PARAM_DESCRIPTION, Default: 10
#' @param wd PARAM_DESCRIPTION, Default: getwd()
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso
#'  \code{\link[e1071]{tune.wrapper}}, \code{\link[e1071]{tune.control}}
#'  \code{\link[randomForest]{randomForest}}, \code{\link[randomForest]{importance}}
#' @rdname sa_rfreg
#' @export
#' @importFrom e1071 tune.randomForest tune.control
#' @importFrom randomForest randomForest importance
sa_rfreg <- function(vpred,
                     custom.name = FALSE,
                     date = "",
                     dataframe.name = "",
                     predictors.name = "",
                     mtry = 1:10,
                     ntree = (1:10)*500,
                     nset = 10, #Number of train/hold sets
                     wd = getwd()
)
{
  # === Output setup ===
  #Output parent folder
  setwd(paste(wd))

  #Output folder
  if(custom.name == FALSE){
    save.name <- readline("Input save.name: RF_df_date_InputVariables_PredictVariable:")
  }
  if (custom.name == TRUE){
    save.name <- paste(dataframe.name, date, predictors.name, vpred, sep = "_")
  }

  dir.create(save.name)
  directory <- paste(wd, save.name, sep = "/")
  print(directory)
  setwd(directory)

  # === Train Random Forest ===
  #RF Loop
  mbest.tab <- data.frame()
  mperf.tab <- data.frame()
  rfimp.tab <- data.frame()
  classmatrix.names <- ""
  print("Initiating RF Code")
  for (i in 1:nset){
    print(i)
    print(Sys.time())

    mbest <- data.frame()
    mperf <- data.frame()
    rfimp <- data.frame()

    tname <- paste("train", i, sep = ".")
    train <- get(tname)
    hname <- paste("hold", i, sep = ".")
    hold <- get(hname)

    print("Initate Tuning")
    #Tuning
    rf = e1071::tune.randomForest(train[,-which(names(train) %in% vpred)], #Remove predicted label from train df
                                  eval(parse(text = paste("train", vpred, sep = "$"))), #Select predicted label (train$vpred)
                                  data = train,
                                  type = "regression", mtry = mtry, ntree = ntree,
                                  tunecontrol = e1071::tune.control(sampling = "cross",cross=5))

    print("Tuning Complete")
    #Model performance table
    name <- rep(tname, length(mtry) * length(ntree))
    mperf <- cbind(name, rf$performances)
    mperf.tab <- rbind(mperf.tab, mperf)





    # === Validate Random Forest ===

    #Test/Hold-Out Performance
    #rfp = random forest with best parameters
    rfp = randomForest::randomForest(train[,-which(names(train) %in% vpred)], #Remove predicted label from train df
                                     eval(parse(text = paste("train", vpred, sep = "$"))), #Select predicted label (train$vpred)
                                     data = train,
                                     type = "regression", mtry = rf$best.parameters$mtry,
                                     ntree = rf$best.parameters$ntree, importance = TRUE)


    rfimp<- randomForest::importance(rfp) %>% data.frame() %>% rownames_to_column(var = "predictors") %>% mutate(set = i)
    rfimp.tab <- rbind(rfimp.tab, rfimp) #rf importance table for export




    #Hold-out set predictions and performance
    pred <- predict(rfp,hold)
    actual <- eval(parse(text = paste0("hold$", vpred)))

    pred.tab <- cbind(pred, actual) %>% data.frame %>% mutate(mpe = ((actual - pred)/actual)*100) %>% mutate(set = i)
    mse <- mean((pred - actual)^2)
    rmse <- sqrt(mse)
    pearson <- cor(pred, actual)
    mape <- mean(((abs(actual - pred))/actual)*100) #mean absolute percent error; better than MPE because +/- absolute values don't cancel out
    mae <- mean((abs(actual-pred))) #mean absolute error

    assign(paste("rf.tab", i, sep = "."), pred.tab) #list of actual/predicted

    #Hold set statistics of best performing model for each train/test set
    name <- tname #Model best parameters & ME/ARI
    mbest <- data.frame(cbind(name, rf$best.parameters$mtry, rf$best.parameters$ntree, mse, rmse, pearson, mape, mae))
    mbest.tab <- data.frame(rbind(mbest.tab,mbest))

    #List of actual/predicted (regression)
    classmatrix.names <- c(classmatrix.names, paste("rf.tab", i, sep = "."))
    print("Hold-Out Assessed")
  } #rf trainings/testing code

  print("RF Complete")

  rf.cm <- c()
  for (i in classmatrix.names){
    print(i)
    if(i != ""){
      cm <- get(i)
      print(cm)
      rf.cm <- rbind(rf.cm, cm)
    }
  } #concatenates actual vs. predicted tables


  #=== Exports ===
  #Save RF optimization data
  colnames(mbest.tab) <- c("set", "mtry", "ntree", "mse", "rmse", "pearson", "mape", "mae")
  write.csv(mbest.tab, "rf_param.csv", row.names = FALSE) #Performance of best performing model for each train/test
  write.csv(mperf.tab, "rf_tuningperformance.csv", row.names = FALSE) #All performances
  write.csv(rf.cm, "rf_predictions.csv", row.names = FALSE) #All predicted vs. actual for all samples
  write.csv(rfimp.tab, "rf_imp.csv", row.names = FALSE) #Variable importance table

  print("CSV Output Complete")

  #Aggregate and save variable importance (average MSE/average Gini)
  varimp2 <- rfimp.tab %>% group_by(predictors) %>%
    summarise_at(vars(X.IncMSE, IncNodePurity), list(name = mean)) %>% arrange(desc(IncNodePurity_name))
  write.csv(varimp2, "imp_aggregate.csv")
  print("Variable importance output complete")

  #Scatterplot of actual vs. predicted values
  tiff(paste0(vpred, "_scatter.tiff"), width = 500, height = 400, res = 100)
  res <- rf.resplot(rf.cm)
  print(res)
  dev.off()

  print("Scatterplot output complete")

  setwd(paste(wd))


  #=== Return results to environment ===
  return(list(varimp = varimp2,
              rf.predictions = rf.cm,
              best.performance = mbest.tab))
}



# Function: Classification RF ---------------------------------------------
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param vpred PARAM_DESCRIPTION
#' @param mtry PARAM_DESCRIPTION, Default: 1:10
#' @param ntree PARAM_DESCRIPTION, Default: (1:10) * 500
#' @param wd PARAM_DESCRIPTION, Default: wd.in
#' @param custom.name PARAM_DESCRIPTION, Default: FALSE
#' @param date PARAM_DESCRIPTION, Default: ''
#' @param dataframe.name PARAM_DESCRIPTION, Default: ''
#' @param predictors.name PARAM_DESCRIPTION, Default: ''
#' @param nset PARAM_DESCRIPTION, Default: 10
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso
#'  \code{\link[e1071]{tune.wrapper}}
#'  \code{\link[randomForest]{importance}}
#' @rdname sa_rf
#' @export
#' @importFrom e1071 tune.randomForest
#' @importFrom randomForest importance
sa_rf <- function(vpred, mtry = 1:10, ntree = (1:10)*500,
                  wd = wd.in,
                  custom.name = FALSE,
                  date = "",
                  dataframe.name = "",
                  predictors.name = "",
                  nset = 10){

  # === Output setup ===
  #Output parent folder
  setwd(paste(wd))

  #Output folder
  if(custom.name == FALSE){
    save.name <- readline("Input save.name: df_RF_date_InputVariables_PredictVariable:")
  }
  if (custom.name == TRUE){
    save.name <- paste(dataframe.name, date, predictors.name, vpred, sep = "_")
  }

  dir.create(save.name)
  directory <- paste(wd, save.name, sep = "/")
  print(directory)
  setwd(directory)


  # === Train Random Forest ===
  #RF loop:
  mbest.tab <- data.frame()
  mperf.tab <- data.frame()
  rfimp.tab <- data.frame()
  classmatrix.names <- ""
  print("Initiating RF Code")
  for (i in 1:nset){
    print(i)
    print(Sys.time())

    mbest <- data.frame()
    mperf <- data.frame()

    tname <- paste("train", i, sep = ".")
    train <- get(tname)
    hname <- paste("hold", i, sep = ".")
    hold <- get(hname)

    print("Initate Tuning")
    #Tuning
    rf = e1071::tune.randomForest(train[,-which(names(train) %in% vpred)], #Remove predicted label from train df
                                  eval(parse(text = paste("train", vpred, sep = "$"))), #Select predicted label (train$vpred)
                                  data = train,
                                  type = "class", mtry = mtry, ntree = ntree,
                                  tunecontrol = tune.control(sampling = "cross",cross=5))

    print("Tuning Complete")
    #Model performance table
    name <- rep(tname, length(mtry) * length(ntree))
    mperf <- cbind(name, rf$performances)
    mperf.tab <- rbind(mperf.tab, mperf)


    # === Validate Random Forest ===
    #rfp = random forest with best parameters
    rfp = randomForest(train[,-which(names(train) %in% vpred)], #Remove predicted label from train df
                       eval(parse(text = paste("train", vpred, sep = "$"))), #Select predicted label (train$vpred)
                       data = train,
                       type = "class", mtry = rf$best.parameters$mtry,
                       ntree = rf$best.parameters$ntree, importance = TRUE)

    pred=predict(rfp,hold,type="class")

    #Model accuracy:
    rf.tab <- table(eval(parse(text = paste("hold", paste(vpred), sep = "$"))),pred)
    diag <- 1-classAgreement(rf.tab)$diag
    crand <- classAgreement(rf.tab)$crand



    assign(paste("rf.tab", i, sep = "."), rf.tab)

    name <- tname #Model best parameters & ME/ARI
    mbest <- data.frame(cbind(name, rf$best.parameters$mtry, rf$best.parameters$ntree, diag, crand))
    mbest.tab <- data.frame(rbind(mbest.tab,mbest))


    classmatrix.names <- c(classmatrix.names, paste("rf.tab", i, sep = "."))

    #Variable importance
    rfimp<- randomForest::importance(rfp) %>% data.frame() %>% rownames_to_column(var = "predictors") %>% mutate(set = i)
    rfimp.tab <- rbind(rfimp.tab, rfimp) #rf importance table for export

    print("Hold-Out Assessed")
  } #rf code
  print("RF Complete")


  #Concatenate confusion matrices:
  rf.cm <- c()
  for (i in classmatrix.names){
    print(i)
    if(i != ""){
      cm <- get(i)
      print(cm)
      rf.cm <- rbind(rf.cm, cm)
    }
  }

  #=== Exports ===
  #Save RF optimization data:
  colnames(mbest.tab) <- c("set", "mtry", "ntree", "ME", "ARI")
  write.csv(mbest.tab, "rf param.csv")
  write.csv(mperf.tab, "rf tuning performance.csv")
  write.csv(rf.cm, "rf cm.csv")
  write.csv(rfimp.tab, "rf_imp.csv", row.names = FALSE) #Variable importance table

  print("CSV Output Complete")



  #OUTPUT: variable importance for the best model/training set
  #rfb = Best random forest model from all X training/test sets (e.g X = 10 sets)
  #Aggregate and save variable importance (average MSE/average Gini)
  varimp2 <- rfimp.tab %P>% group_by(predictors) %P>% summarise_at(vars(MeanDecreaseAccuracy, MeanDecreaseGini), mean) %>% arrange(desc(MeanDecreaseGini))
  print(head(varimp2))
  write.csv(varimp2, "imp_aggregate.csv")
  print("Aggregated variable importance output complete")



  print("Variable importance complete")

  setwd(paste(wd))

  #=== Return results to environment ===
  return(list(varimp = varimp2,
              rf.cm = rf.cm,
              best.performance = mbest.tab))
}




# Functions: Variable importance formatting --------------------------------
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param df PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname varimport_import
#' @export
varimport_import <- function(df)
  #edit varimp df imported from CSV
{
  df.out <- df %>% rename(predictors = X)
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param rf.type PARAM_DESCRIPTION
#' @param metric PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname select.varimp
#' @export
select.varimp <- function(rf.type, metric)
  #Select variable importance metric (gini, mse, mda)
{
  #Unit test: rf.type and metric inputs compatible.
  try(if((rf.type == "class" & metric == "mse") | (rf.type == "reg" & metric == "mda"))
    stop("rf.type and metric incompatabile. Ensure you are using the correct variable importance metric for your random forest type.")
  )



  if(rf.type == "reg"){
    if(metric == "mse"){
      y <- "X.IncMSE_name"
      ylab <- "Increased Mean Square Error"
    }

    if(metric == "gini"){
      y <- "IncNodePurity_name"
      ylab <- "Increased node purity"
    }
  }

  if(rf.type == "class"){
    if(metric == "mda"){
      y <- "MeanDecreaseAccuracy"
      ylab <- "Decreased Accuracy"
    }

    if(metric == "gini"){
      y <- "MeanDecreaseGini"
      ylab <- "Increased node purity"
    }

  }
  return(c(y, ylab))
}


# Function: Variable Importance -------------------------------------------
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param varimp PARAM_DESCRIPTION
#' @param rf.type PARAM_DESCRIPTION
#' @param metric PARAM_DESCRIPTION, Default: 'gini'
#' @param selection_type PARAM_DESCRIPTION, Default: 'random_top'
#' @param top PARAM_DESCRIPTION, Default: 20
#' @param xlab PARAM_DESCRIPTION, Default: NA
#' @param extract.names.df PARAM_DESCRIPTION, Default: NA
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso
#'  \code{\link[dplyr]{count}}
#' @rdname varimport_plot
#' @export
#' @importFrom dplyr count
varimport_plot <- function(varimp, #dataframe of variable importance
                           rf.type, #One of: "class" or "reg"
                           metric = "gini", #One of: "mse", "gini", or "mda)
                           selection_type = "random_top", #One of: random_top, top
                           top = 20, #number of "top" predictor variables to plot
                           xlab = NA, #Labels to use for plotting. NA will use R column names. Accepts a character vector or "extract" instructions.
                           extract.names.df = NA  #provide df if xlab = "extract". 1st col is desired x label names. 2nd col is how columns appear in df input.
)
  #Lollipop plot of variable importance

  #Varimp df structure:
  #Predictor variable column must be named ("predictors")
  #Use varimp import


  #Selection_type:
  #top = top 20 (or specified) predictor variable (PV)
  #random_top = #Take top 20 (or specified) predictor variables (PV), bottom 5 (specified*0.25) PV, and 15 (specified*0.75) additional random PV for visualization

{
  #=== X and y variables ===
  select <- select.varimp(rf.type, metric)
  y <- select[1]
  ylab <- select[2]

  print(y) #Check y
  print(ylab) #Check y label

  x <- "predictors" #do not change




  #=== Variable importance dataframe ===
  if(selection_type == "random_top")
  {
    toppred <- varimp %>% slice_max(order_by = get(y), n = top)
    minpred <- varimp %>% slice_min(order_by = get(y), n = round(top*0.25, digits = 0))
    randompred <- varimp %>% filter(!(predictors %in% toppred$predictors | predictors %in% minpred$predictors)) %>%
      slice_sample(n = round(top*0.75, digits = 0))

    varimp <- varimp %>% filter(predictors %in% toppred$predictors | predictors %in% minpred$predictors | predictors %in% randompred$predictors)
  }

  if(selection_type == "top")
    #Take top 20 taxa (or specified)
  {
    varimp <- varimp %>% slice_max(order_by = get(y), n = top)
  }

  varimp <- varimp %>% arrange(desc(y)) #re-order variables





  #=== Predictor variable labels ===
  xlab.edit <- ifelse(is.na(xlab),
                      FALSE, #No user-input x labels provided. Use varimp column names.
                      TRUE) #edit predictor labels according to user-input


  if(!(is.na(xlab)) & xlab == "extract")
    #extract method
  {
    #Unit test: extract.names.df provided
    if(class(extract.names.df) == "logical")
    {
      if(extract.names.df[1] %>% is.na){
        print("Error: Provide dataframe to extract variables.")
        print("Running UNDECLARED to kill script.")
        UNDECLARED()
      }
    }

    #Align extract.names.df to varimp order
    xlab.extract <- varimp %>% rename(column_names = "predictors") %>%
      left_join(extract.names.df, by = "column_names")

    #Unit test: Top_pred and extract.names.df labels align
    qc <- (xlab.extract$column_names == varimp$predictors) %>% data.frame() %>% rename(qc = ".") %>% dplyr::count(qc) %>% unique
    print(qc)
    if(qc[1] == FALSE){
      print("Error: Extract names dataframe labels do not align with predictor variable order")
      print("Running UNDECLARED to kill script.")
      UNDECLARED()
    }


    xlab <- xlab.extract
  }

  #Unit test: are there enough xlabels supplied?
  if(is.na(xlab) == FALSE){
    if(nrow(xlab) != nrow(varimp))
    {
      print("Error: Mismatched number of x labels and predictor variables")
      print("Running UNDECLARED to kill script.")
      UNDECLARED()
    }
  }



  #=== Plot ===
  p <- ggdotchart(varimp, x = paste(x), y = paste(y),
                  sorting = "ascending",
                  ggtheme = theme_pubr(),
                  add = "segment",
                  ylab = ylab,
                  xlab = "") +
    {if(xlab.edit == TRUE)
      scale_x_discrete("", breaks = xlab$column_names, labels = xlab$xlabel_names)
    } +
    coord_flip() +
    theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5))

  return(p)
}




# Functions: Density plot loop  --------------------------------------------
#Density plot function:
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param df PARAM_DESCRIPTION
#' @param x PARAM_DESCRIPTION
#' @param dvar PARAM_DESCRIPTION
#' @param tsize PARAM_DESCRIPTION, Default: 24
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname gdens
#' @export
gdens <- function(df, #name of dataframe
                  x, #region to plot
                  dvar, #group by (e.g genotype)
                  tsize = 24
)
{
  ggplot(df, aes_string(x = x, fill = dvar)) +
    geom_density(alpha=0.5) +
    theme_minimal_hgrid(12) +
    theme(axis.text.y=element_blank(),
          axis.title.y=element_blank(),
          axis.ticks.y=element_blank(),
          text = element_text(size=tsize))
}





#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param df PARAM_DESCRIPTION
#' @param response PARAM_DESCRIPTION
#' @param tsize PARAM_DESCRIPTION, Default: 24
#' @param top PARAM_DESCRIPTION, Default: 20
#' @param varimp PARAM_DESCRIPTION
#' @param rf.type PARAM_DESCRIPTION
#' @param metric PARAM_DESCRIPTION
#' @param scale PARAM_DESCRIPTION, Default: 0.8
#' @param ncol PARAM_DESCRIPTION, Default: 2
#' @param xlab PARAM_DESCRIPTION, Default: NA
#' @param extract.names.df PARAM_DESCRIPTION, Default: NA
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso
#'  \code{\link[cowplot]{plot_grid}}
#' @rdname multi.density.plot
#' @export
#' @importFrom cowplot plot_grid
multi.density.plot <- function(df, #dataframe containing predictor variables and response variable by subject
                               response, #response variable (text)
                               tsize = 24,
                               top = 20, #number of "top" predictor variables to plot
                               varimp, #variable importance dataframe imported using
                               rf.type, #One of: "class" or "reg"
                               metric, #One of: "mse", "gini", or "mda
                               scale = 0.8, #cowplot scale graphs
                               ncol = 2, #cowplot number of cols
                               xlab = NA, #one of: NA, "extract", "blank", or list of custom x labels listed according to plot order
                               extract.names.df = NA #provide df if xlab = "extract". 1st col is desired x label names. 2nd col is how columns appear in df input.
)

#extract.names.df formatting:
#1st column = xlabel_names
#X labels to appear in density plot
#2nd column = column_names
#How the variable appears in R columns

{
  #Set variable importance metric
  select <- select.varimp(rf.type, metric)
  metric.name <- select[1]

  top_pred <- varimp %>% slice_max(order_by = get(metric.name), n = top) %>% pull(predictors)

  # === Set x labels ===
  #No x labels
  if(xlab == "blank" & !(is.na(xlab))){
    xlab <- rep("", length(top_pred))
  }

  #X labels extraccted from df
  if(xlab == "extract" & !(is.na(xlab)))
    #extract method
  {
    #Unit test: extract.names.df provided
    if(!(is.na(extract.names.df)))
    {
      if(extract.names.df[1] %>% is.na){
        print("Error: Provide dataframe to extract variables.")
        print("Running UNDECLARED to kill script.")
        UNDECLARED()
      }
    }

    #Align extract.names.df to top_pred order
    xlab.extract <- top_pred %>% data.frame %>% rename(column_names = ".") %>%
      left_join(extract.names.df, by = "column_names")

    #Unit test: Top_pred and extract.names.df labels align
    qc <- (xlab.extract$column_names == top_pred) %>% data.frame() %>% rename(qc = ".") %>% count(qc) %>% unique
    if(qc[1] == FALSE){
      print("Error: Extract names dataframe labels do not align with predictor variable order")
      print("Running UNDECLARED to kill script.")
      UNDECLARED()
    }


    #Extract labels
    xlab <- xlab.extract %>% pull(xlabel_names)

  }

  #Unit test: are there enough xlabels supplied?
  if(is.na(xlab) == FALSE){
    if(length(xlab) != length(top_pred)){
      print("Error: Mismatched number of x labels and predictor variables")
      print("Running UNDECLARED to kill script.")
      UNDECLARED()
    }
  }

  #Multiple plotting:
  dens.df <- data.frame()
  plotlist <- list()
  for (i in 1:length(top_pred))
  {
    j <- top_pred[i]
    print(i)
    print(j)

    p <- gdens(df = df,
               x = j,
               dvar = response,
               tsize = tsize) + theme(legend.position = "none") +
      {if(is.na(xlab) == FALSE)
        xlab(xlab[i])}


    assign(paste("dens", i, sep ="."), p)
    plotname <- paste("dens", i, sep =".")
    dens.df <- rbind(dens.df, plotname)
    plotlist[[i]] <- p
  }

  #Extract legend
  plotlist[[length(plotlist) + 1]] <- gdens(df = df, x = top_pred[1], dvar = response, tsize = tsize) %>% get_legend()

  #Plot grid
  p_grid <- cowplot::plot_grid(plotlist = plotlist, scale = scale, ncol = ncol)
  return(p_grid)
}



# Classification Matrix Aggregation ---------------------------------------
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param cm PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname cm.aggregate
#' @export
cm.aggregate <- function(cm #cm object output from random forest functions
)
  #First column of input cm must be named "X"
  #This is the default behaviour upon importing cm output from random forest functions
  #User must not change column names from original output
{
  cm.groups <- cm %>% select(-X) %>% colnames
  cm.out <- cm %>% group_by(X) %>% summarise(across(cm.groups, ~ sum(.)))
  return(cm.out)
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param agg.cm PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname cm.prop
#' @export
cm.prop <- function(agg.cm #aggregated confusion matrix. Output from cm.aggregate function.
)
  #Expresses confusion matrix as the percentage of each class being classified by each other class
{
  cm.groups <- agg.cm %>% select(-X) %>% colnames
  cm.out <- agg.cm %>% mutate(across(cm.groups, ~ (./sum(., na.rm = TRUE)) %>% round(digits = 2)))
  return(cm.out)
}








# Mean and Stdev of Model Performance Indicators  -------------------------
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param rf.type PARAM_DESCRIPTION
#' @param performance.table PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname mean_var.rfmodel
#' @export
mean_var.rfmodel <- function(rf.type, #"class" or "reg"
                             performance.table #mbest.tab output
){
  if(rf.type == "class"){
    mp <- performance.table %>% mutate(across(.cols = c("ME", "ARI"), .fns = as.numeric)) %>%
      summarise(across(.cols = c("ME", "ARI"),
                       .fns = list(Mean = mean, SD = sd)
      ))
    model.performance.out <- data.frame(mean = c(mp$ME_Mean, mp$ARI_Mean),
                                        sd = c(mp$ME_SD, mp$ARI_SD))
    row.names(model.performance.out) <- c("ME", "ARI")
  }

  if(rf.type == "reg"){
    mp <- performance.table %>% mutate(across(.cols = c("mse", "rmse", "pearson", "mape", "mae"), .fns = as.numeric)) %>%
      summarise(across(.cols = c("mse", "rmse", "pearson", "mape", "mae"),
                       .fns = list(mean = mean, SD = sd)))
    print(mp)
    model.performance.out <- data.frame(mean = c(mp$mse_mean, mp$rmse_mean, mp$pearson_mean, mp$mape_mean, mp$mae_mean),
                                        sd = c(mp$mse_SD, mp$rmse_SD, mp$pearson_SD, mp$mape_SD, mp$mae_SD)
    )

    row.names(model.performance.out) <- c("mse", "rmse", "pearson", "mape", "mae")
  }
  return(model.performance.out)

}



