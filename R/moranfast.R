# Modifications to function from moranfast package, authored by Matthew Cooper

#' Calculate Moran's I statistic quickly
#'
#' @param x a numeric vector of values with a spatial distribution
#' @param lat a numeric vector indicating the latitude for each element of x
#' @param lon a numeric vector indicating the longitude for each element of x
#' @param alternative a string specifying the alternative hypothesis to be tested against, must be one of "two.sided", "less", or "greater"
#'
#' @description
#' Calculates Moran's I statistic, measuring spatial autocorrelation, for a
#' vector of values distributed across space, with locations indicated by
#' latitude/longitude. Utilizes improvements originally developed by Matt Cooper
#' to optimize the calculation for big data applications. It is memory optimized
#' as it calculates the distance matrix on-the-fly as needed for the I statistic
#' calculation, rather than defining the entire spatial weights matrix first, then
#' calculating the global Moran's I. Calculations are performed in C++ using Rcpp
#' to speed up the process over implementations calculating Moran's I that are
#' purely written in R.
#'
#' @export
moranfast <- function(x, lat, lon, alternative='two.sided'){
  res <- calc_moran(x, lat, lon)

  names(res) <- c('observed', 'expected', 'sd')
  res <- as.list(res)

  alternative <- match.arg(alternative, c("two.sided", "less", "greater"))
  pv <- pnorm(res$observed, mean = res$expected, sd = res$sd)
  if (alternative == "two.sided"){
    if (res$observed <= -1/(length(x) - 1)){
      pv <- 2 * pv
    }else{
      pv <- 2 * (1 - pv)
    }
  }
  if (alternative == "greater"){
    pv <- 1 - pv
  }

  res[['p.value']] <- pv

  return(res)
}
