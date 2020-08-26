## ---------------------------------------------------------
##                NON-EXPORTED HELPERS
## ---------------------------------------------------------
# Used in query
update_marginal_tables <- function(am, z) {
  lapply(am, function(x) {
    if (!is.null(attr(x, "vars"))) {
      # if (setequal(attr(x, "vars"), c("V50", "V58"))) browser()
      idx <- match(attr(x, "vars"), colnames(z))
      el  <- paste(z[idx], collapse = "")
      if (el %in% names(x)) {
        x[el] <- x[el] + 1L
      } # Else do nothing since a 0 and a 1 evaluates to the same
    } else { # Handling empty separators
      x <- x + 1L
    }
    return(x)
  })
}

# used in plot.multiple_models
extract_model_simulations <- function(models) {
  if (!inherits(models, "multiple_models")) stop("`models` needs to be an object returned from `fit_multiple_models`")
  sims <- lapply(seq_along(models), function(m) {
    data.frame(Deviance = models[[m]]$sims,
      response = names(models)[m],
      stringsAsFactors = FALSE)
  }) 
  do.call(rbind, sims)
}

# used in plot.multiple_models
make_observation_info <- function(models) {
  zdevs <- sapply(models, function(m) m$dev)
  zpvs  <- sapply(models, function(m) m$pv)
  data.frame(devs = zdevs, pvals = zpvs, response = names(models))
}


## ---------------------------------------------------------
##                   EXPORTED HELPERS
## ---------------------------------------------------------

#' Print outlier model
#'
#' A print method for \code{outlier_model} objects
#'
#' @param x A \code{outlier_model} object
#' @param ... Not used (for S3 compatability)
#' @export
print.outlier_model <- function(x, ...) {
  cls <- paste0("<", paste0(class(x)[-length(class(x))], collapse = ", "), ">")
  cat(
    "\n --------------------------------",
    "\n  Simulations:",         length(x$sims),
    "\n  Variables:",           ncol(x$A),
    "\n  Observations:",        nrow(x$A),
    "\n  Estimated mean:",      round(x$mu, 2),
    "\n  Estimated variance:",  round(x$sigma, 2),
    paste0("\n  ", cls),
    "\n --------------------------------\n"
  )
}

print.outlier_model_approx <- function(x, ...) {
  print.outlier_model(x$m, ...)
}


#' Print outlier
#'
#' A print method for \code{outlier} objects
#'
#' @param x A \code{outlier} object
#' @param ... Not used (for S3 compatability)
#' @export
print.outlier <- function(x, ...) {
  cls <- paste0("<", paste0(class(x)[-length(class(x))], collapse = ", "), ">")
  cat(
    "\n --------------------------------",
    "\n  Simulations:",         length(x$sims),
    "\n  Variables:",           ncol(x$A),
    "\n  Observations:",        nrow(x$A),
    "\n  Estimated mean:",      round(x$mu, 2),
    "\n  Estimated variance:",  round(x$sigma, 2),
    "\n    ---------------------------  ",
    "\n  Critical value:", x$cv,
    "\n  Deviance:", x$dev,
    "\n  P-value:", x$pval,
    "\n  Alpha:", x$alpha,
    paste0("\n  ", cls),
    "\n --------------------------------\n"
  )
}


#' Plot Distribution of Simulated Deviances
#'
#' A plot method to show the approximated deviance distribution
#' @param x An object returned from \code{fit_outlier}
#' @param sig_col Color of the significance level area (default is red)
#' @param ... Extra arguments; see details.
#' @details The dotted line represents the observed deviance of the observation under the hypothesis
#' and the colored (red is default) area under the graph represents the significance level.
#' Thus, if the dotted line is to the left of the colored area, the hypothesis that the observation
#' is an outlier cannot be rejected. Notice however, if there is no dotted line, this simply means,
#' that the observed deviance is larger than all values and it would disturb the plot if included.
#'
#' No extra arguments \code{...} are implement at the moment.
#' @export
plot.outlier <- function(x, sig_col = "#FF0000A0", ...) {
  # args <- list(...)
  # Old base approach:
  # graphics::hist(x$sims, breaks = 30, xlab = "Deviance",  freq = FALSE, main = " ")
  # dat <- data.frame(Deviance = x$sims, y = "")
  dat <- with(stats::density(x$sims), data.frame(x, y))
  dat$.region <- x$cv
  dat$.dev    <- x$dev
  p <- ggplot2::ggplot(data = dat, mapping = ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_line() +
    ggplot2::geom_ribbon(data = subset(dat, x >= .region),
      ggplot2::aes(ymax = y),
      ymin   = 0,
      fill   = sig_col,
      colour = sig_col,
      alpha  = 0.7) + 
    ggplot2::geom_ribbon(data = subset(dat, x <= .region),
      ggplot2::aes(ymax = y),
      ymin   = 0,
      fill   = "#A0A0A0A0",
      colour = "#A0A0A0A0",
      alpha  = 0.7)
  p <- p + ggplot2::theme_bw() + ggplot2::ylab("") + ggplot2::xlab("Deviance")
  if (dat$.dev[1] < max(x$sims)) {
    p <- p + ggplot2::geom_vline(ggplot2::aes(xintercept = .dev), linetype = "dotted")
  }
  return(p)
}

#' Plot Deviance of Multiple Models
#'
#' A plot method to show the approximated deviance distribution of multiple models
#' @param x A \code{multiple_models} object returned from a called to \code{fit_multiple_models}
#' @param sig_col Color of the significance level area (default is red)
#' @param ... Extra arguments. See details.
#' @details The dotted line represents the observed deviance of the observation under the hypothesis
#' and the colored (red is default) area under the graph represents the significance level.
#' Thus, if the dotted line is to the left of the colored area, the hypothesis that the observation
#' is an outlier cannot be rejected. Notice however, if there is no dotted line, this simply means,
#' that the observed deviance is larger than all values and it would disturb the plot if included.
#'
#' No extra arguments \code{...} are implement at the moment.
#' @export
plot.multiple_models <- function(x, sig_col = "#FF0000A0", ...) {
  z_dev_pval <- make_observation_info(x)
  dat        <- extract_model_simulations(x)
  p <- ggplot2::ggplot(dat, ggplot2::aes(x = Deviance, y = response))
  p <- p + ggridges::stat_density_ridges(ggplot2::aes(fill=factor(..quantile..)),
    geom      = "density_ridges_gradient",
    calc_ecdf = TRUE,
    quantiles = c(1 - x[[1]]$alpha, 1)
  )
  p <- p + ggplot2::scale_y_discrete(limits = names(x))
  p <- p + ggplot2::scale_fill_manual(
    name  = "Significance level",
    values = c("#A0A0A0A0", sig_col, sig_col),
    labels = c("", "(ms[[1]]$alpha, 1]", "")
  )
  
  for (i in 1:nrow(z_dev_pval)) {
    max_dev_i <- max(dat$Deviance)
    dev <- z_dev_pval[i, "devs"]
    if (dev < max_dev_i) { # Dont plot the dotted line if it is larger than all deviances
      linet  <- "dotted"
      df_seg <- data.frame(x1 = dev, x2 = dev, y1 = i , y2 = i + 1)
      p <- p + ggplot2::geom_segment(ggplot2::aes(x = x1,
        y    = y1,
        xend = x2,
        yend = y2
      ),
      linetype = linet,
      size     = 1,
      color    = "black",
      data     = df_seg
      )
    }
  }
  p <- p + ggplot2::theme_bw() +
    ggplot2::ylab("") +
    ggplot2::xlab("Deviance") +
    ggplot2::theme(legend.position = "none")
  return(p)
}


#' Plot the probability mass function of the deviances
#'
#' A plot method to show the pmf of the approximated pmf of \code{D(Y)}
#'
#' @param x A \code{outlier_model} object
#' @param alpha Significance level
#' @param sig_col Color of the significance level area (default is red)
#' @param ... Not used (for S3 compatibility)
#' @export
plot.outlier_model <- function(x, alpha = 0.05, sig_col = "#FF0000A0", ...) {
  # graphics::hist(x$sims, breaks = 30, xlab = "Deviance",  freq = FALSE, main = " ")
  dat <- with(stats::density(x$sims), data.frame(x, y))
  dat$.region <- critval(x, alpha)
  p <- ggplot2::ggplot(data = dat, mapping = ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_line() +
    ggplot2::geom_ribbon(data = subset(dat, x >= .region),
      ggplot2::aes(ymax = y),
      ymin   = 0,
      fill   = sig_col,
      colour = sig_col,
      alpha  = 0.7) + 
    ggplot2::geom_ribbon(data = subset(dat, x <= .region),
      ggplot2::aes(ymax = y),
      ymin   = 0,
      fill   = "#A0A0A0A0",
      colour = "#A0A0A0A0",
      alpha  = 0.7)
  p <- p + ggplot2::theme_bw() + ggplot2::ylab("") + ggplot2::xlab("Deviance")
  return(p)
}

#' Query Outlier Tests
#'
#' Given an outlier model, an outlier test for each row in the new data is conducted
#'
#' @param m A \code{outlier_model} object
#' @param z data.frame with new data
#' @param alpha The significance level
#' @details Use \code{fit_outlier_model} to obtain \code{m}
#' @export
query <- function(m, z, alpha = 0.05) UseMethod("query")

#' @rdname query
#' @export
query.outlier_model <- function(m, z, alpha = 0.05) {
  if (!neq_empt_chr(disc(m$vs))) stop("Not implemented for 'pure_cont' models yet.")
  if (!m$approx) stop("m is not an approximated model. Use approx = TRUE in fit_outlier_model")
  z_disc     <- z[, disc(m$vs)]
  new_models <- lapply(1:nrow(z), function(i) {
    zi   <- z_disc[i, , drop = FALSE]
    cms  <- update_marginal_tables(m$cms, zi)
    sms  <- update_marginal_tables(m$sms, zi)
    sims_disc <- 2*(.map_dbl(m$sims_disc_raw, TY, cms, sms) - Hx_(nrow(m$A)+1L))
    sims <- m$sims_mixed + sims_disc
    # Update the model
    mi       <- m
    mi$A     <- rbind(mi$A, z)
    mi$sims  <- sims
    mi$sims_disc_raw <- NULL
    mi$mu    <- mean(sims)
    mi$sigma <- stats::var(sims)
    mi$cdf   <- stats::ecdf(sims)
    mi$cms   <- cms
    mi$sms   <- sms
    mi$pa    <- m$pa
    mi$vs    <- m$vs
    return(mi)
  })
  lapply(1:nrow(z), function(i) {
    mi   <- new_models[[i]]
    cv   <- critval(mi, alpha)
    devi <- deviance(mi, z[i, ])
    pi   <- pval(mi, devi)
    new_outlier(mi, devi, pi, cv, alpha)
  })
}

#' Query Outlier Tests - version 2
#'
#' Given an outlier model, an outlier test for each row in the new data is conducted
#'
#' @inheritParams query
#' @details Use \code{fit_outlier_model} to obtain \code{m}
#' @export
query2 <- function(m, z, alpha = 0.05) UseMethod("query2")

#' @rdname query2
#' @export
query2.outlier_model <- function(m, z, alpha = 0.05) {
  cv <- critval(m, alpha)
  lapply(1:nrow(z), function(i) {
    mi   <- m
    # mi$A <- rbind(mi$A, z[i, ])  
    devi <- deviance(mi, z[i, ])
    new_outlier(mi, devi, pval(mi, devi), cv, alpha)
  })
}

#' Empirical distribution function
#'
#' The empirical cdf of \code{D(Y)}
#'
#' @param x A \code{outlier_model} object
#' @param ... Not used (for S3 compatibility)
#' @export
cdf <- function(x, ...) UseMethod("cdf")

#' @rdname cdf
#' @export
cdf.outlier_model <- function(x, ...) return(x$cdf)

#' P-value
#'
#' Calculate the p-value for obtaining \code{ty_new} under \code{H_0}
#'
#' @param x A \code{outlier_model} object
#' @param dz The deviance of the observation \code{z}.
#' @param ... Not used (for S3 compatibility)
#' @details The value \code{dz} can be obtained used the \code{deviance} function.
#' @seealso \code{\link{deviance}}
#' @export
pval <- function(x, dz, ...) UseMethod("pval")

#' @rdname pval
#' @export
pval.outlier_model <- function(x, dz, ...) return(1 - x$cdf(dz))


#' Critical value
#'
#' Calculate the critical value for test statistic under \code{H_0}
#'
#' @param m A \code{outlier_model} object
#' @param alpha Significance level (between \code{0} and \code{1})
#' @details The value \code{dz} can be obtained used the \code{deviance} function.
#' @seealso \code{\link{deviance}}
#' @export
critval <- function(m, alpha = 0.05) UseMethod("critval")

#' @rdname critval
#' @export
critval.outlier_model <- function(m, alpha = 0.05) {
  stats::uniroot(function(x) pval(m, x) - alpha,
    interval = range(m$sims),
    extendInt = "yes",
    tol = 0.0001)$root
}

#' Mean
#'
#' Estimated mean of deviance statistic \code{T(Y)}
#'
#' @param x A \code{outlier_model} object
#' @param ... Not used (for S3 compatibility)
#' @export
mean.outlier_model <- function(x, ...) return(x$mu)

#' Variance
#'
#' Estimated variance of the deviance statistic \code{T(Y)}
#'
#' @param x A \code{outlier_model} object
#' @param ... Not used (for S3 compatibility)
#' @export
variance <- function(x) UseMethod("variance")

#' @rdname variance
#' @export
variance.outlier_model <- function(x, ...) return(x$sigma)
