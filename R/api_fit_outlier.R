#' Outlier model
#'
#' A model based on decomposable graphical models for outlier detection
#'
#' @param A data,frame with the new observation of interest appended
#' @param adj Adjacency list or gengraph object of a decomposable graph without the
#' observation of interest. See details
#' @param nsim Number of simulations
#' @param ncores Number of cores to use in parallelization
#' @param validate Logical. See details.
#' @param approx Logical. The user should rarely set this manually.
#' @param hom Logical indicating wheter the model is homogeneous or not if \code{A} has mixed variables
#' @details It is assumed that all cell values in \code{A}, for all variables,
#' are represented as a single character. If \code{validate} is \code{TRUE} this is checked.
#' If cell values are not single characters, one may exploit the \code{to_single_chars} function.
#'
#' The \code{adj} object is most typically found using \code{fit_graph}. But the user can supply
#' an ajacency list, just a named \code{list}, of their own choice if needed.
#'
#' The parameter \code{approx} is used to output the correct type of discrete simulations
#' whether or not an approximation is to be used. \code{fit_outlier} is an exact test
#' and calls \code{fit_outlier_model} with \code{approx = TRUE}.
#' @seealso \code{\link{fit_graph}}, \code{\link{query}}, \code{\link{fit_outlier}}
#' @export
fit_outlier_model <- function(A,
                              adj,
                              nsim       = 10000,
                              ncores     = 1,
                              validate   = TRUE,
                              approx     = TRUE,
                              hom        = FALSE
                              ) {
  env <- environment()
  if (inherits(adj, "gengraph")) adj <- adj_lst(adj)
  M    <- nrow(A)
  vs   <- verts(A)
  type <- type_of_verts(vs)

  Ad   <- if (neq_empt_chr(disc(vs))) as.matrix(A[, disc(vs), drop = FALSE]) else NULL # A

  if (validate && type != "pure_cont") {
    if (!only_single_chars(Ad)) {
      stop("All discrete values in A must be represented as a single character. See '?to_single_chars'") 
    }
  }

  sims_disc_raw <- NULL
  sims_mixed    <- 0L
  mp <- marginals_and_parents(adj, type, vs, Ad)

  sims <- if (type == "pure_cont") {
    .sim_internal_pure_cont(mp$pa, M, nsim)
  } else if (type == "pure_disc") {
    if (approx) {
      sim_disc <- .sim_internal_pure_disc(Ad, mp$cms, mp$sms, type = "pure_disc_raw", nsim, ncores)
      assign("sims_disc_raw", sim_disc, envir = env)
      2 * (.map_dbl(sim_disc, TY, mp$cms, mp$sms) - Hx_(M))
    } else {
      .sim_internal_pure_disc(Ad, mp$cms, mp$sms, type = "pure_disc", nsim, ncores)
    }
  } else {
    sim_disc   <- .sim_internal_pure_disc(Ad, mp$cms, mp$sms, type = "pure_disc_raw", nsim, ncores)
    sim_mixed  <- .sim_internal_mixed(mp$pa, sim_disc, hom) # parallelize ?
    if (approx) {
      assign("sims_disc_raw", sim_disc, envir = env)
      assign("sims_mixed",    sim_mixed, envir = env)
    }
    sim_disc  <- 2 * (.map_dbl(sim_disc, TY, mp$cms, mp$sms) - Hx_(M))
    sim_mixed + sim_disc
  }
  mu    <- mean(sims)
  sigma <- stats::var(sims)
  cdf   <- stats::ecdf(sims)
  m <- new_outlier_model(
    A,
    sims,
    sims_mixed,
    sims_disc_raw,
    mu,
    sigma,
    cdf,
    mp$cms,
    mp$sms,
    mp$pa,
    type,
    vs,
    hom,
    approx
  )
  return(m)
}

#' Fit Outlier
#'
#' Outlier test in decomposable graphical models
#'
#' @param A Data without the new observation \code{z} appended (data.frame)
#' @param z New observation with same colnames as A (data.frame)
#' @param adj Adjacency list or gengraph object of a decomposable graph
#' without the observation of interest. See details
#' @param alpha The significance level
#' @param nsim Number of simulations
#' @param ncores Number of cores to use in parallelization
#' @param validate Logical. If true, it checks if \code{A} has only single character values and converts it if not.
#' @param hom Logical indicating whether the model is homogeneous or not if \code{A} has mixed variables
#' @details It is very important to note, that it is assumed that \code{z} is
#' not appended to \code{A}. It may be, that \code{z} actually belongs to \code{A} and we want
#' to test if it is an outlier. But then one should pass the subdata \code{A[-zi, ]} where
#' \code{zi} is the row index of \code{z} in \code{A}. The reason for this is, that under the
#' null hypothesis, \code{z} is believed to stem from \code{A} and so in the procedure
#' \code{z} is appended. But if \code{z} is already in \code{A} this may lead to a slightly
#' different result since \code{z} was in reality appended twice.
#'
#' The \code{adj} object is most typically found using \code{fit_graph}. But the user can supply
#' an adjacency list, just a named \code{list}, of their own choice if needed.
#' @seealso \code{\link{fit_outlier_model}}, \code{\link{fit_graph}}
#' @export
fit_outlier <- function(A,
                        z,
                        adj,
                        alpha    = 0.05,
                        nsim     = 10000,
                        ncores   = 1,
                        validate = TRUE,
                        hom      = FALSE) {

  # TODO:
  # - let z be null if z is not a new outlier just like in the molic package

  if (!identical_colnames(A, z)) stop("Variables in A and z is not in agreement!")
  if (inherits(adj, "gengraph")) adj <- adj_lst(adj)

  Az   <- rbind(A, z)
  vs   <- verts(Az)
  type <- type_of_verts(vs)

  ## Azd  <- if (neq_empt_chr(disc(vs))) Az[, disc(vs), drop = FALSE] else NULL
  if (validate && neq_empt_chr(disc(vs))) {
    if (!only_single_chars(Az[ ,disc(vs), drop = FALSE])) {
      stop("All discrete values in A must be represented as a single character. See '?to_single_chars'") 
    }
  }

  m     <- fit_outlier_model(Az, adj, nsim, ncores, validate = FALSE, approx = FALSE, hom)
  dev_z <- deviance(m, z)
  m     <- new_outlier(m, dev_z, pval(m, dev_z), critval(m, alpha), alpha)
  return(m)
}
