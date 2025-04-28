# Modifications to ME function from spatialreg package
# To support limited set of eigenvectors
# ME function from spatial reg originally copyright 2005-2008
# by Roger Bivand and Pedro Peres-Neto (from Matlab)

#' Performs Moran eigenvector spatial filtering
#'
#' @param formula a symbolic description of the regression model to be fit
#' @param data a data frame containing the variables in the model
#' @param lat latitude column in data
#' @param lon longitude column in data
#' @param family the error distribution and link function to be used in the regression model
#' @param weights an optional vector of weights to be used in the model fitting process
#' @param offset an optional a priori known component to be included in the linear predictor during fitting
#' @param na.action a function to specify behavior for missing data (default options("na.action")), can also be na.omit or na.exclude
#' @param listw a spatial weights list, for example created by nblistw
#' @param alpha a stopping rule - no further eigenvectors will be added to model once p-value for residual autocorrelation exceeds alpha
#' @param verbose if TRUE prints messages indicating progress in model fitting
#' @param zero.policy if FALSE stop with error for any empty neighbor sets, if TRUE permit the weights list to be formed with zero-length weights vectors
#' @param n_eigs limit testing of which moran eigenvectors to include in regression model to the n_eigs eigenvectors of the spatial weights matrix with the largest magnitude
#'
#' @description
#' Quickly performs Moran eigenvector spatial filtering (MESF). MESF involves
#' removing/reducing spatial autocorrelation present in the residuals of
#' generalized linear models (GLMs). MESF searches eigenvectors of a transformed
#' spatial weights matrix (a doubly centered matrix). One-by-one, eigenvectors
#' are added to the right-hand side of the GLM and Moran's I statistic (a measure
#' of spatial autocorrelation) is calculated for the residuals of the newly fitted
#' model. The eigenvector which produces the lowest Moran's I is then permanently
#' added to the model and the search process repeated to choose additional
#' eigenvectors to add to the model. This process continues until the p-value
#' for Moran's I exceeds the specified stopping threshold alpha.
#'
#' This implementation optimizes MESF for bigger data applications in several ways:
#' 1. Calculates Moran's I in a faster and more memory efficient way using the
#' moranfast() function implemented in this package. This function calculates
#' Moran's I on the fly rather than by calculating and storing in memory a new
#' distance matrix for every new set of residuals produced, and uses Rcpp/C++
#' to speed up calculation.
#' 2. Creating the transformed spatial weights matrix (the eigenvectors
#' of which will be considered for inclusion in the model) using faster matrix
#' computations from Rfast package which runs in parallel in C++.
#' 3. Reducing the set of eigenvectors over which to conduct a brute force search
#' to only the n_eigs largest eigenvectors of the transformed spatial weights
#' matrix. This speeds up the process by only computing the smaller number of
#' eigenvectors we are going to consider (using the RSpectra package) and by
#' requiring us to test fewer eigenvectors in the model re-fitting stage.
#'
#' @export
mesf <- function(formula, data=list(), lat, lon, family = gaussian,
                 weights, offset, na.action=na.fail,
                 listw=NULL, alpha=0.05, verbose=NULL,
                 zero.policy=NULL, n_eigs = 100) {

  #### running in parallel
  cores <- parallel::detectCores()
  if (is.null(cores)) {
    parallel <- "no"
  } else {
    parallel <- ifelse (spatialreg::get.mcOption(), "multicore", "snow")
  }

  ncpus <- ifelse(is.null(cores), 1L, cores-1)
  cl <- NULL
  if (parallel == "snow") {
    cl <- spatialreg::get.ClusterOption()
    if (is.null(cl)) {
      parallel <- "no"
      warning("no cluster in ClusterOption, parallel set to no")
    }
  }
  par_boot_args <- list(parallel=parallel, ncpus=ncpus, cl=cl)
  #####

  ##### handle missing arguments & defaults
  if (is.null(zero.policy))
    zero.policy <- get("zeroPolicy", envir = spatialreg:::.spatialregOptions)
  stopifnot(is.logical(zero.policy))

  if (is.null(verbose)) verbose <- get("verbose", envir = spatialreg:::.spatialregOptions)
  stopifnot(is.logical(verbose))

  # argument handling copied from stats:::glm
  #	call <- match.call()
  if (is.character(family))
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family)) family <- family()
  if (is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }
  #####

  #### extract X & Y
  #     	if (missing(data)) data <- environment(formula)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "weights", "offset", "na.action"),
             names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())

  mt <- attr(mf, "terms")
  Y <- stats::model.extract(mf, "response") # extract X and Y
  X <- stats::model.matrix(mt, mf)
  #####

  #### handle more specialized args - weights, offset, listw
  weights <- stats::model.weights(mf)
  if (!is.null(weights) && any(weights < 0))
    stop("negative weights not allowed")

  offset <- stats::model.offset(mf)
  if (!is.null(offset) && length(offset) != NROW(Y))
    stop("number of offsets should equal number of observations")

  na.act <- attr(mf, "na.action")
  if (is.null(listw)) stop("listw required")
  if (!is.null(na.act)) {
    subset <- !(1:length(listw$neighbours) %in% na.act)
    listw <- subset(listw, subset, zero.policy=zero.policy)
  }
  #####



  #### fit initial regression model with no eigenvector filtering
  if (verbose) cat("fitting initial model \n")
  glm_fit <- stats::glm.fit(x=X, y=Y, weights=weights, offset=offset,
                     family=family)
  glm_res <- glm_fit$y - glm_fit$fitted.values
  if (verbose) cat("found initial model \n")
  #####


  ### calculate initial moran statistic
  # fast calculation of moran statistic using rcpp
  if (verbose) cat("getting initial moran i \n")
  mRes <- moranfast(glm_res, lat, lon)
  pIZ <- mRes$p.value
  tres <- c(NA, mRes$observed, pIZ)
  if (pIZ > alpha) stop("base correlation larger than alpha")
  if (verbose) cat("Initial Moran I: ", mRes$observed, " Initial p-value: ",pIZ,"\n", sep="")
  ###

  #### create doubly centered symmetric distance matrix
  # and calculate n_eigs largest eigenvectors from this matrix
  if (verbose) cat("converting to full matrix form \n")
  Wmat <- spdep::listw2mat(listw) # convert to full matrix form
  n <- ncol(Wmat)

  if (verbose) cat("doubly centering CWC matrix \n")
  R <- matrix(0,ncol=n, nrow=n) + Rfast::rowmeans(Wmat)
  if (verbose) cat("got rowmeans, now calculating col means \n")
  C <- matrix(0,ncol=n, nrow=n) + Rfast::colmeans(Wmat, parallel=TRUE, cores=ncpus)
  if (verbose) cat("got colmeans, now taking transpose \n")
  C <- Rfast::transpose(C)
  if (verbose) cat("getting mean of whole matrix \n")
  rowsums <- Rfast::rowsums(Wmat, parallel=TRUE, cores=ncpus)
  grandsum <- sum(rowsums)
  if (verbose) cat("combining row and colmean info to add to grand mean \n")
  CWC <- Wmat - R - C + grandsum/n
  rm(Wmat, R, C, rowsums, grandsum) # clean up local environ to reduce memory usage

  if(!isSymmetric(CWC)){
    if (verbose) cat("making CWC matrix symmetric \n")
    CWC <- 0.5*(CWC + Rfast::transpose(CWC))
  }

  if (verbose) cat("got CWC matrix, getting eigs \n")
  eV <- RSpectra::eigs_sym(CWC, k=n_eigs, which="LM")$vectors
  if (verbose) cat("finished calculating eigs \n")
  rm(CWC)
  ####

  ### fit model with each eigenvector to find the first eigenvector to add
  if (verbose) cat("looking for best eigenvector for eig1 \n")
  iZ <- numeric(n_eigs)
  for (i in 1:n_eigs) {
    iX <- cbind(X, eV[,i])
    i_glm <- stats::glm.fit(x=iX, y=Y, weights=weights, offset=offset,
                     family=family)
    glm_res <- i_glm$y - i_glm$fitted.values
    iZ[i] <- moranfast(glm_res, lat, lon)$observed
  }
  min_iZ <- which.min(abs(iZ))
  ####

  ### add first eigenvector to list of xs and refit model
  X <- cbind(X, eV[, min_iZ])
  glm_fit <- stats::glm.fit(x=X, y=Y, weights=weights, offset=offset,
                     family=family)
  glm_res <- glm_fit$y - glm_fit$fitted.values

  mRes <- moranfast(glm_res, lat, lon)
  pIZ <- mRes$p.value
  used <- rep(FALSE, n_eigs)
  used[min_iZ] <- TRUE
  min_v <- min_iZ
  if (verbose) cat("eV[,", min_iZ, "], I: ", mRes$observed, " Expected: ",
                   mRes$expected, ", pr(ZI): ", pIZ, "\n", sep="")
  tres <- rbind(tres, c(min_iZ, mRes$observed, pIZ))
  #####

  #### look for all the rest of the eigs to use
  while (pIZ <= alpha) {
    if (all(used)==TRUE){
      if (verbose) cat("all eigs used \n")
      break
    }
    for (i in 1:n_eigs) {
      if (used[i]==FALSE) {
        iX <- cbind(X, eV[,i])
        i_glm <- stats::glm.fit(x=iX, y=Y, weights=weights,
                         offset=offset, family=family)
        glm_res <- i_glm$y - i_glm$fitted.values

        iZ[i] <- moranfast(glm_res, lat, lon)$observed
      } else iZ[i] <- NA
    }
    min_iZ <- which.min(abs(iZ))
    X <- cbind(X, eV[, min_iZ])
    glm_fit <- stats::glm.fit(x=X, y=Y, weights=weights, offset=offset,
                       family=family)
    glm_res <- glm_fit$y - glm_fit$fitted.values

    mRes <- moranfast(glm_res, lat, lon)
    pIZ <- mRes$p.value
    used[min_iZ] <- TRUE
    min_v <- c(min_v, min_iZ)
    if (verbose) cat("eV[,", min_iZ, "], I: ", mRes$observed, " Expected: ",
                     mRes$expected, ", pr(ZI): ", pIZ, "\n", sep="")
    tres <- rbind(tres, c(min_iZ, mRes$observed, pIZ))
  }

  if(verbose) cat("got all eigs \n")
  sv <- eV[,min_v, drop=FALSE]
  colnames(sv) <- paste("vec", min_v, sep="")
  colnames(tres) <- c("Eigenvector", "ZI", "pr(ZI)")
  rownames(tres) <- 0:(nrow(tres)-1)
  res <- list(selection=tres, vectors=sv)
  if (!is.null(na.act))
    res$na.action <- na.act
  class(res) <- "Me_res"
  res
}

print.Me_res <- function(x, ...) {
  print(x$selection)
}

fitted.Me_res <- function(object, ...) {
  if (is.null(object$na.action)) {
    res <- object$vectors
  } else {
    omitted_rows <- unname(object$na.action)
    res <- matrix(as.numeric(NA), ncol=ncol(object$vectors),
                  nrow=length(omitted_rows)+nrow(object$vectors))
    i <- j <- k <- 1
    while (i <= nrow(res)) {
      if (j <= length(omitted_rows) && i == omitted_rows[j]) {
        i <- i+1
        j <- j+1
      } else {
        res[i,] <- object$vectors[k,]
        i <- i+1
        k <- k+1
      }
    }
  }
  res
}
