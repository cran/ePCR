% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ext.R
\name{conforminput}
\alias{conforminput}
\title{Conform the dimensions of a new input data matrix to a readily fitted PEP or PSP object}
\usage{
conforminput(object, newx)
}
\arguments{
\item{object}{A readily fitted PSP or PEP object}

\item{newx}{A data matrix or a data.frame which to expand}
}
\value{
An expanded data matrix for which the dimensions conform to the regression coefficients in the PSP or PEP
}
\description{
Conform the dimensions of a new input data matrix to a readily fitted PEP or PSP object
}
\author{
Teemu Daniel Laajala \email{teelaa@utu.fi}
}
