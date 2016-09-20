
#' Read MAVEN Data
#'
#' Reads a MAVEN datafile from directory
#'
#' @param directory a pathway pointing to an input file which is either .xlsx or tab-delimited txt/tsv
#'
#' @return a list:
#' metabolite_info: information specific to an individual metabolite
#' metabolomics_data: metabolite abundances (ion counts) in each sample in tidy format
#'
#' @export
#'
#' @examples
#' directory = system.file("extdata", "maven_axon.xlsx", package = "fluxr")
#' read_MAVEN_data(directory)
read_MAVEN_data <- function(directory, ...){

  possible_id_variables <- c("label", "metaGroupId", "groupId", "goodPeakCount", "medMz", "medRt", "maxQuality", "note",
                             "compound", "compoundId", "expectedRtDiff", "ppmDiff", "parent")

  if(!file.exists(directory)){
    stop(directory, ": file does not exist")
  }

  if(!grepl('(.xlsx)|(.txt)|(.tsv)$', directory)){
    stop("Only .xlsx or tab-delimited text (.txt / .tsv) currently supported")
  }

  if(grepl('.xlsx$', directory)){
    maven_file <- readxl::read_excel(directory)
  }else{
    maven_file <- readr::read_tsv(directory)
  }

  metabolite_info <- maven_file %>%
    dplyr::select_(.dots = colnames(maven_file)[colnames(maven_file) %in% possible_id_variables]) %>%
    dplyr::mutate(metabolite_row = 1:n()) %>%
    dplyr::select(metabolite_row, dplyr::everything())

  metabolomics_data <- maven_file %>%
    dplyr::select_(.dots = colnames(maven_file)[!(colnames(maven_file) %in% possible_id_variables)]) %>%
    dplyr::mutate(metabolite_row = 1:n()) %>%
    tidyr::gather("sample", "IC", -metabolite_row) %>%
    # convert ion counts of zero to NA
    dplyr::mutate(IC = ifelse(IC == 0, NA, IC))

  maven_data <- list()
  maven_data$metabolite_info <- metabolite_info
  maven_data$metabolomics_data <- metabolomics_data
  maven_data

  }

#' Normalize Metabolomics
#'
#' Using a robust metabolomics normalization rescale ion counts based on the consensus signal of each sample.
#'
#' This normalization was reported in Kamphorst et al. 2015, see \url{https://github.com/shackett/Pancreatic_tumor_metabolomics}.
#'
#' @param metabolomics_data a tidy data.frame generated from \code{\link{read_MAVEN_data}}
#'
#' @return metabolomics_data with an added column containing normalized ion counts
#'
#' @export
#'
#' @examples
#' directory = system.file("extdata", "maven_axon.xlsx", package = "fluxr")
#' normalize_metabolomics(read_MAVEN_data(directory)$metabolomics_data)
normalize_metabolomics <- function(metabolomics_data){

  metabolomics_data %>%
    dplyr::filter(!is.na(IC)) %>%
    # find the median of each metabolite across samples
    dplyr::group_by(metabolite_row) %>%
    dplyr::mutate(median_abund = median(IC)) %>%
    # calculate scaling factors relative to median metabolite
    dplyr::ungroup() %>%
    dplyr::mutate(relative_met_conc = IC/median_abund) %>%
    # find how much more abundant a sample is on average across metabolite
    dplyr::group_by(sample) %>%
    dplyr::mutate(scaling_factor = median(relative_met_conc)) %>%
    # rescale ion counts based on scaling_factor
    dplyr::mutate(IC_normalized = IC / scaling_factor) %>%
    # return original columns and IC_normalized
    dplyr::select_(.dots = c(colnames(metabolomics_data), "IC_normalized"))

  }


