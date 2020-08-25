## ---------------------------------------------------------
##                NON-EXPORTED HELPERS
## ---------------------------------------------------------
parents_formula <- function(y, pc) {
  if (neq_empt_chr(pc)) {
    return(stats::as.formula(paste(y, paste(pc, collapse = " + "), sep = "~")))
  } else {
    return(stats::as.formula(paste(y, " ~ 1L")))
  }
}

reduce_data_to_parent_cells <- function(d, zpd) {
  # zpd : data.rame with nrow = 1. Column names are the discrete parents of a continuous variable in the perfect ordering
  dzpd <- d[, colnames(zpd), drop = FALSE]
  d[apply(dzpd, 1, identical, unlist(zpd)), ]
}

.sum_eta <- function(cell_data, u, v) {
  # cell_data: list of data.frames
  # u,v: characters (not constrained to singletons)

  u_singleton <- length(u) == 1L
  v_singleton <- length(v) == 1L  
  
  se <- lapply(cell_data, function(x) {

    scaled_data <- scale(x[, unique(c(u, v)), drop = FALSE], scale = FALSE)

    out <- if (u_singleton && v_singleton) {
      scaled_data[, u, drop = TRUE] * scaled_data[, v, drop = TRUE]
    } else {
      lapply(1:nrow(scaled_data), function(i) {
        outer(scaled_data[i, u], scaled_data[i, v])
      })
    }

    Reduce("+", out)
  })
  return(se)
}


.SSE <- function(A, y, pa_d, pa_c) {
  # A: All data for variables y, pa_d and pa_c;
  # y: The continous dependent variable
  # pa_d: Names of the discrete parents (character(0) if none)
  # pa_c: Names of the continuous parents (character(0) if none)
  if (neq_empt_chr(pa_d)) {
    A$lvls <- apply(A[, pa_d, drop = FALSE], 1, paste, collapse = "")
    cell_data <- split(A, A$lvls)
    if (neq_empt_chr(pa_c)) { # Both discrete and cont. parents
      Scc <- Reduce("+", .sum_eta(cell_data, pa_c, pa_c))
      Scy <- Reduce("+", .sum_eta(cell_data, pa_c, y))
      Syc <- t(Scy)
      Syy <- Reduce("+", .sum_eta(cell_data, y, y))  ## FIX!!! OUTER IS OVERKILL HERE!
      SSE <- Syy - Syc %*% solve(Scc) %*% Scy
      return(as.vector(SSE))
    } else { # Only discrete parents
      Syy <- Reduce("+", .sum_eta(cell_data, y, y))
      return(as.vector(Syy))
    }
  } else {
    Al  <- list(A) # .sum_eta assumes the first argument to be a list of data.frames
    if (neq_empt_chr(pa_c)) { # Only continous parents
      Scc <- Reduce("+", .sum_eta(Al, pa_c, pa_c))
      Scy <- Reduce("+", .sum_eta(Al, pa_c, y))
      Syc <- t(Scy)
      Syy <- Reduce("+", .sum_eta(Al, y, y))
      SSE <- Syy - Syc %*% solve(Scc) %*% Scy
      return(as.vector(SSE))
    } else { # No parents
      Syy <- Reduce("+", .sum_eta(Al, y, y))
      return(as.vector(Syy))
    }
  }
}

## ---------------------------------------------------------
##                   EXPORTED HELPERS
## ---------------------------------------------------------

#' Calculate deviance
#'
#' This function calculates the affine value \code{-2 log} likelihood-ratio statistic; also called the deviance statistic
#'
#' @param x A \code{outlier_model} object
#' @param z An observation (data.frame or named character vector)
#' @param ... Not used (for S3 compatability)
#' @export
deviance <- function(x, z, ...) {
  UseMethod("deviance")
}

#' @rdname deviance
#' @export
deviance.pure_disc <- function(x, z,...) {
  # FIX: if y is not a character an error is produced that stops the R process
  # - it also needs the same colnames as x$A etc... fix!
  2 * (TY(unlist(z), x$cms, x$sms) - Hx_(nrow(x$A))) # D(y)
}

#' @rdname deviance
#' @export
deviance.pure_cont <- function(x, z, ...) {
  # z is/should already be appended to x$A in outlier_model
  # if (!inherits(x, "pure_cont")) stop("outlier object not of type 'pure_cont'")
  # z_in_A    <- row_match(z, x$A)
  # if (is.na(z_in_A)) stop("z was never included in the data when fit_outlier was called!")
  A       <- x$A[-nrow(x$A), ]# x$A[-z_in_A, ]
  M       <- nrow(A) + 1L
  pa_c    <- attr(x$pa, "pc")
  cont    <- names(pa_c)
  n_paj_c <- length(pa_c)
  -M * sum(.map_dbl(1:n_paj_c, function(j) {
    paj_c <- pa_c[[j]]
    yj    <- cont[j]
    pf    <- parents_formula(yj, paj_c)
    # This can, in principle go wrong, if k > n_i0 (lm will raise an error if so!)
    r0 <- summary(stats::lm(pf, rbind(A, z)))$residuals
    ra <- summary(stats::lm(pf, A))$residuals
    Qj <- sum(ra^2) / sum(r0^2)
    log(Qj)
  }))
}


#' @rdname deviance
#' @export
deviance.mixed <- function(x, z, ...) { # make a "check" argument for z_in_A
  rownames(z) <- NULL
  if (!identical_colnames(x$A, z)) stop("z does not agree with data!")
  if (inherits(z, "tbl_df"))     z <- as.data.frame(z)
  if (inherits(x$A, "tbl_df")) x$A <- as.data.frame(x$A)
  if (x$hom) {
    return(.dev_mixed_hom(x, z, ...))
  } else {
    return(.dev_mixed_inhom(x, z, ...))
  }
}

.dev_mixed_hom <- function(x, z, ...) {
  pa_d    <- attr(x$pa, "pd")
  pa_c    <- attr(x$pa, "pc")
  pa_n    <- attr(x$pa, "nd")
  cont    <- names(pa_c)
  # Assume that z is the last element for approx = FALSE!
  A <- if (x$approx) x$A else x$A[-nrow(x$A), ]
  M       <- nrow(A) # Assuming that z is included
  dev_Q <- sum(.map_dbl(seq_along(pa_c), function(j) {
    paj_d  <- pa_d[[j]]
    paj_c  <- pa_c[[j]]
    zpaj_d <- z[, paj_d, drop = FALSE]
    yj     <- cont[j]
    Aj     <- A[, c(yj, paj_d, paj_c)]
    SSEj0  <- .SSE(Aj, yj, paj_d, paj_c)
    SSEj   <- .SSE(Aj[-nrow(Aj), ], yj, paj_d, paj_c) ## Assume z is the last row in A!!
    Qj     <- log(SSEj / SSEj0)
    return(Qj)
  }))
  dev_QD <- 2 * (TY(unlist(z[ , x$vs$disc]), x$cms, x$sms) - Hx_(M)) # D(y)
  return(dev_Q + dev_QD)
}

.dev_mixed_inhom <- function(x, z, ...) {
  ## z is/should already be appended to x$A in fit_outlier or in query
  # Assume that z is the last element for approx = FALSE!
  A       <- if (x$approx) x$A else x$A[-nrow(x$A), ]
  pa_d    <- attr(x$pa, "pd")
  pa_c    <- attr(x$pa, "pc")
  pa_n    <- attr(x$pa, "nd")
  cont    <- names(pa_c)
  n_paj_c <- length(pa_c)
  dev_Q  <- sum(.map_dbl(1:n_paj_c, function(j) {
    paj_n  <- pa_n[[j]]
    paj_d  <- pa_d[[j]]
    paj_c  <- pa_c[[j]]
    yj     <- cont[j]
    zpaj_d <- z[, paj_d, drop = FALSE]
    n_zpaj <- if (neq_empt_chr(paj_d)) paj_n[paste(unlist(zpaj_d), collapse = "")] else paj_n
    Aj     <- if (neq_empt_chr(paj_d)) reduce_data_to_parent_cells(A, zpaj_d) else A
    pf     <- parents_formula(yj, paj_c)
    ## ---------------------------------------------------------
    ##               STUDENTIZED VERSION
    ## ---------------------------------------------------------
    ## tictoc::tic()
    ## ma      <- stats::lm(pf, Aj)
    ## mu_hat  <- predict(ma, newdata = z)
    ## yj_0    <- z[, yj, drop = TRUE]
    ## ypaj_0  <- if (neq_empt_chr(paj_c)) as.numeric(z[, paj_c]) else numeric(0L)
    ## xx      <- c(1, ypaj_0)
    ## xx_vcov <- 1 + sum(t(xx) %*% vcov(ma) %*% xx)
    ## se_res  <- summary(ma)$sigma * sqrt(xx_vcov)
    ## res     <- yj_0 - mu_hat
    ## dfj     <- n_zpaj - length(paj_c) - 1
    ## F_      <- (res / se_res)^2
    ## Qj      <- 1 - (F_/dfj) / (1 + F_/dfj)
    ## tictoc::toc()
    ## ---------------------------------------------------------

    # Find the correlation matrix and pick the two regressors
    # with the highest correlation. Then remove the one, that
    # correlates the least with the response variable. Keep doing
    # this until Residual standard error !is.nan

    ## cor(Aj[, c(yj, paj_c)])[-1, -1]
    ## cor(Aj[, c(yj, paj_c)])[1, ][-1]
    
    r0 <- summary(stats::lm(pf, rbind(Aj, z)))$residuals
    ra <- summary(stats::lm(pf, Aj))$residuals
    Qj <- sum(ra^2) / sum(r0^2)
    # if (is.nan(log(Qj))) browser()
    -n_zpaj * log(Qj)
  }))
  dev_QD <- 2 * (TY(unlist(z[ , x$vs$disc]), x$cms, x$sms) - Hx_(nrow(A))) # D(y)  
  return(dev_Q + dev_QD)
}

## ---------------------------------------------------------
##                 TEST HEART DATA
## ---------------------------------------------------------
## library(dplyr)

## d <- heart_orig_wide_without_VF %>%
##   select(-case)

## g <- fit_mixed_graph(d, TRUE)
## plot(gengraph(d, "gen", g), vertex.size = 5)


## tests <- lapply(1:nrow(d), function(i) {
##   print(i)
##   mi <- fit_outlier(d[-i, ], d[i, ], g, validate = FALSE, ncores = 3, nsim = 2000) 
##   return(mi)
## })

## heart_orig_wide_without_VF[sapply(tests, function(x) x$pval <= 0.05), ] %>% select(case)

## ---------------------------------------------------------
##                    EXAMPLE 1
## ---------------------------------------------------------
## library(dplyr)

## g = list(
##   A = c("B", "X", "Y"),
##   B = c("A", "Y"),
##   X = c("A", "Y"),
##   Y = c("A", "X", "B")
## )

## .lvls <- list(
##   A = c("0", "1"),
##   B = c("0", "1")
## )

## d <- cgr_sim(g, .lvls, ns = 1500,
##   alpha_rate = 10,
##   beta_rate  = 10,
##   sigma_rate = 10,
##   cell_rate  = 5) %>%
##   as_tibble()

## m <- fit_outlier_model(d, g, nsim = 1000, validate = FALSE, approx = TRUE)
## deviance(m, d[1,])


## ## ## ---------------------------------------------------------
## ## ##              COMPLETE GRAPH
## ## ## ---------------------------------------------------------

## g <- make_complete_graph(letters[1:10])
## .lvls <- structure(replicate(4, c("0", "1", "3", "4", "5"), FALSE), names = names(g)[1:4])

 
## d <- cgr_sim(g, .lvls, ns = 500,
##   alpha_rate = 10,
##   beta_rate  = 10,
##   sigma_rate = 10,
##   cell_rate  = 1) %>%
##   as_tibble()

## m <- fit_outlier_model(d, g, nsim = 10000, validate = FALSE, approx = TRUE)
## plot(m)
## deviance(m, d[1,])
