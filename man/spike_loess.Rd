% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RAnanlysis.R
\name{spike_loess}
\alias{spike_loess}
\title{Normalize data using loess}
\usage{
spike_loess(data, subset)
}
\arguments{
\item{data}{A numeric matrix of raw counts.}

\item{subset}{A vector of indices to subset the data for normalization.}
}
\value{
A loess-normalized numeric matrix.
}
\description{
This function performs loess normalization on the dataset.
}
\examples{
# Example usage
normalized_data <- spike_loess(data, subset = c(1:10))
}
