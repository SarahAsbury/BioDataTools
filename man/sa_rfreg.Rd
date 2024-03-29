% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RandomForest_Functions.R
\name{sa_rfreg}
\alias{sa_rfreg}
\title{sa_rfreg}
\usage{
sa_rfreg(
  vpred,
  custom.name = FALSE,
  date = "",
  dataframe.name = "",
  predictors.name = "",
  mtry = 1:10,
  ntree = (1:10) * 500,
  nset = 10,
  wd = getwd()
)
}
\arguments{
\item{vpred}{Predictor variable as string.}

\item{custom.name}{(Optional) \code{TRUE} or \code{FALSE}. Will the user input a custom name for the export directory? If FALSE, user must ineractively provide one in the console. \strong{Recommend to set \code{custom.fun  = TRUE} TRUE for automation}. Default: FALSE}

\item{date}{(Optional). If \code{custome.name = FALSE} provide date in the input in a format suitable for folder/directory name.}

\item{dataframe.name}{(Optional). If \code{custom.name = FALSE} provide short name of the dataset.}

\item{predictors.name}{(Optional). If \code{custom.name = FALSE} provide short name of the dataset.}

\item{mtry}{(Optional) Set range for optimization of mtry in randomForest function (predictor variables to use at each split). \strong{Recommend to determine using mtry.guide()} Default: 1:10}

\item{ntree}{(Optional) Set range for optimization of ntrees in randomForest function, Default: (1:10)*500}

\item{nset}{Number of train/hold sets to try.}

\item{wd}{}
}
\value{

}
\description{
Run random forest regression for train/hold sets labeled train.1 to train.i and hold.1 to hold.i. This function is built to run with Split dataset function in the rf_standard wrapper. If the user does not provide a custom directory name (requires input of date, dataframe.name, and predictors.name), this function requires interaction in the R console.
}
\seealso{
[randomForest()], [mtry.guide()], [split.ratio()],[rf_standard()]
}
