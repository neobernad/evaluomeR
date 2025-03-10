% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dataPreprocessing.R
\name{determineNumberOfFactors}
\alias{determineNumberOfFactors}
\title{Determine Number of Factors}
\usage{
determineNumberOfFactors(dataset)
}
\arguments{
\item{dataset}{A data frame from which factors are to be extracted.
The function will ignore the "Description" column if it exists.}
}
\value{
An integer representing the number of factors to retain based on the
        scree test results.
}
\description{
This function determines the optimal number of factors to extract from a dataset
using the eigenvalue-based scree test and parallel analysis.
}
