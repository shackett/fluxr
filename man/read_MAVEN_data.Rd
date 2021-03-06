% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/metabolomics.R
\name{read_MAVEN_data}
\alias{read_MAVEN_data}
\title{Read MAVEN Data}
\usage{
read_MAVEN_data(directory, ...)
}
\arguments{
\item{directory}{a pathway pointing to an input file which is either .xlsx or tab-delimited txt/tsv}
}
\value{
a list:
metabolite_info: information specific to an individual metabolite
metabolomics_data: metabolite abundances (ion counts) in each sample in tidy format
}
\description{
Reads a MAVEN datafile from directory.
}
\examples{
# load MAVEN data from an excel file (tab-delimited .txt and .tsv are also accepted)
directory = system.file("extdata", "maven_axon.xlsx", package = "fluxr")
read_MAVEN_data(directory)

# loading data with isotopologue abundances
directory = system.file("extdata", "maven_axon_isotopes.xlsx", package = "fluxr")
read_MAVEN_data(directory)
}

