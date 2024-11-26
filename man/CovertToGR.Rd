% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RAnanlysis.R
\name{CovertToGR}
\alias{CovertToGR}
\title{Convert Data Frame to Genomic Ranges}
\usage{
CovertToGR(data, ID = "Geneid")
}
\arguments{
\item{data}{A data frame containing genomic data.}

\item{ID}{A character string specifying the column name for the ID field (default: "Geneid").}
}
\value{
A GRanges object representing the genomic data.
}
\description{
This function converts a data frame with genomic coordinates into a GRanges object.
}
\examples{
# Example usage
gr <- CovertToGR(data, ID = "Geneid")
}
