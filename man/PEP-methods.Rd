% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pep.R
\name{PEP-methods}
\alias{PEP-methods}
\alias{print,PEP-method}
\alias{predict,PEP-method}
\title{PEP-methods}
\usage{
\S4method{print}{PEP}(x, ...)

\S4method{predict}{PEP}(object, type = "response", newx, x.expand)
}
\arguments{
\item{x}{Generic x}

\item{...}{Additional custom parameters passed on}

\item{object}{PEP-ensemble model object}

\item{type}{Type of prediction; either "response" or "ensemble"}

\item{newx}{New data matrix}

\item{x.expand}{A function that may expand (i.e. extract features) from the input data matrix. By default this will be the default x.expand saved in the S4-slot of the first ensemble member. If the user wishes to omit this functionality, setting this parameter to 'x.expand = as.matrix' does not expand the input data matrix. Notice that if the user has manually called the 'conforminput' function for the newx-data, it is no longer necessary to expand the data matrix here.}
}
\description{
PEP-methods

print.PEP: Print a PEP object to the console

predict.PEP: Predict for a novel patient from current PEP-ensemble
}
