#' @name columbus_sp
#' @title Columbus crime data by neighborhood
#'
#' @description Data on crime rates, household income, and home values for
#' different neighborhoods in Columbus, OH from 1980. Modified from columubs data
#' available in the sp package.
#'
#' @format This data is a SpatialPointsDataFrame with 49 observations of 6
#' variables:
#'    \itemize{
#'      \item{id: ID number of a particular neighborhood}
#'      \item{crime: residential burglaries and vehicle thefts per thousand
#'              households in the neighborhood}
#'      \item{home_value: housing value in 1,000 USD}
#'      \item{income: household income in 1,000 USD}
#'      \item{open_space: amount of open space in neighborhood}
#'      \item{no_plumbing: percentage housing units without plumbing}
#'      \item{dist_to_cbd: distance to central buisness district}
#'      \item{lat: approximate latitude for center of neighborhood}
#'      \item{lon: approximate longitude for center of neighborhood}
#'    }
"columbus_sp"
