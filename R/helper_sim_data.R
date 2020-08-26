## ---------------------------------------------------------
##                 NON-EXPORTED HELPERS
## ---------------------------------------------------------
parents <- function(g) {
  # out: the parents in a perfect ordering
  # g  : a pure discrete graph
  ps <- mcs_marked(g, disc_vars = names(g))$ps # need the marked version to get "ps"
  structure(
    lapply(seq_along(ps), function(i) {
      setdiff(ps[[i]], names(ps)[i])
    }),
    names = names(ps)
  )
}

p_array <- function(pas, lvls, cell_rate = 0.5) {
  # pas: parents of a perfect ordering 
  p_arr <- lapply(seq_along(pas), function(i) {
    p         <- names(pas)[i]
    pasi      <- pas[[i]]
    dim_names <- if (neq_empt_chr(pasi)) c(pasi, p) else p
    dim_vec   <- .map_int(lvls[dim_names], length)
    pi_arr    <- cbind(expand.grid(lvls[c(pasi, p)], stringsAsFactors = FALSE), prob = 0L)

    for (x in pasi) { # x in character(0) convieniently skips the loop
      pi_arr <- pi_arr[order(pi_arr[, x]), ]
    }

    reps_of_p_lvls <- 1:(nrow(pi_arr) / dim_vec[p])
    pi_arr[, "prob"] <- as.vector(
      sapply(reps_of_p_lvls, function(x) {
        px <- stats::rbeta(unname(dim_vec[p]), 2, cell_rate) # 2 is magic!
        px / sum(px)
      })
    )
    
    return(pi_arr)
  })
  structure(p_arr, names = names(pas))
}

sim_pure_disc_data <- function(g, lvls, nsim = 1000, cell_rate = 0.5) {
  # out: Data from a discrete markov random field
  # ----------------------------------
  # g   : a pure discrete graph
  # lvls: levels of the discrete variables - a named list
  pas  <- parents(g)
  parr <- p_array(pas, lvls, cell_rate)
  y    <- matrix("", nrow = nsim, ncol = length(names(g)), dimnames = list(NULL, names(g)))
  ## PARALLELIZE!!!??
  for (i in 1:nsim) {
    yi <- y[i, ]
    for (p in parr) {
      curr_var <- colnames(p)[(ncol(p)-1)]
      if (ncol(p) == 2L) {
        yi[curr_var] <- sample(p[, curr_var], 1L, prob = p[, "prob"])
      } else {
        curr_pas <- pas[[curr_var]]
        curr_val <- yi[curr_pas]
        for (cp in curr_pas) {
          p <- p[p[, cp, drop = TRUE] == curr_val[cp], ]
        }
        yi[curr_var] <- sample(p[, curr_var], 1L, prob = p[, "prob"]) 
      }
    }
    y[i, ] <- yi
  }
  return(y)
}

coefs_pure_cont <- function(pa,
                   alpha_rate = 1,
                   beta_rate  = 1,
                   sigma_rate = 5) {
  n_pa  <- length(pa)
  list(
    alpha = structure(stats::runif(n_pa, -alpha_rate, alpha_rate), names = names(pa)),  
    beta  = structure(stats::runif(n_pa, -beta_rate,  beta_rate),  names = names(pa)),  
    sigma = structure(stats::runif(n_pa, sigma_rate, 5*sigma_rate), names = names(pa)) 
  )
}


sim_pure_cont_data <- function(pa, ns, coefs) {
  empty_obs <- structure(vector("numeric", length(pa)), names = names(pa))
  alpha <- coefs[["alpha"]]
  beta  <- coefs[["beta"]]
  sigma <- coefs[["sigma"]]
  Ac <- lapply(1:ns, function(i) {
    y <- empty_obs
    for (k in names(y)) {
      pa_c <- attr(pa, "pc")[[k]]
      a <- alpha[k]
      b <- beta[k]
      b <- if (neq_empt_chr(pa_c)) seq_along(pa_c) * b else 0L # betas are just linear related then
      s <- sigma[k]
      y[k] <- a + sum(b * y %[chr% pa_c) + stats::rnorm(1, 0, s)
    }
    y
  })
  as.data.frame(do.call(rbind, Ac), stringsAsFactors = FALSE)
}

coefs_mixed <- function(nd,
                   hom,
                   alpha_rate = 1,
                   beta_rate  = 1,
                   sigma_rate = 5,
                   beta_hom   = 1L,
                   sigma_hom  = 5L) {
  M     <- sum(nd[[1]])
  alpha <- lapply(nd, function(x) x / M * stats::runif(length(x), -alpha_rate, alpha_rate))
  beta  <- if (hom) {
    structure(replicate(length(nd), list(beta_hom)), names = names(nd))
  } else {
    lapply(nd, function(x) {
      x / M * stats::runif(length(x), -beta_rate, beta_rate)
    })
  }
  sigma <- if (hom) {
    structure(replicate(length(nd), list(sigma_hom)), names = names(nd))
  } else {
    lapply(nd, function(x) {
      x / M * stats::runif(length(x), 1, sigma_rate) #  # 1 is magic!
    })      
  }
  list(alpha = alpha, beta = beta, sigma = sigma)
}

sim_mixed_data <- function(Ad, pa, ns, coefs, hom) {
  empty_obs <- structure(vector("numeric", length(pa)), names = names(pa))
  alpha <- coefs[["alpha"]]
  beta  <- coefs[["beta"]]
  sigma <- coefs[["sigma"]]
  Ac <- lapply(1:ns, function(i) {
    z <- Ad[i, ]
    y <- empty_obs
    for (k in names(y)) {
      pa_d  <- attr(pa, "pd")[[k]]
      pa_c  <- attr(pa, "pc")[[k]]
      n_pad <- attr(pa, "nd")[[k]]
      z_pa  <- paste(z[1L, pa_d, drop = TRUE], collapse = "")
      a <- alpha[[k]][z_pa]
      b <- if (hom) beta[[k]]  else beta[[k]][z_pa]
      b <- if (neq_empt_chr(pa_c)) seq_along(pa_c) * b else 0L # betas are just linear related then
      s <- if (hom) sigma[[k]] else sigma[[k]][z_pa]
      y[k] <- a + sum(b * y %[chr% pa_c) + stats::rnorm(1, 0, s)
    }
    y
  })
  as.data.frame(do.call(rbind, Ac), stringsAsFactors = FALSE)
}


#' Simulate Observations From a Conditional Gaussian Regression
#'
#' @param g Adjacency list of a decomposable (mixed) graph
#' @param lvls Levels of the discrete variables. A named list.
#' @param hom Logical indicating if the model is homogeneous or not
#' @param ns Number of simulations
#' @param nc Number of cores to simulate the discrete part
#' @param alpha_rate Control for intercept on cell level
#' @param beta_rate Control for slopes on cell level
#' @param sigma_rate Control for variances on cell level
#' @param cell_rate Control discrete cell probabilities
#' @param beta_hom Control for fixed slope if \code{hom = TRUE}
#' @param sigma_hom Control for fixed variance if \code{hom = TRUE}
#' @details cell_rate, beta_hom and sigma_hom are not used in the pure_cont case
#' and when lvls is NULL it is assumed that g is of type pure_cont
#' @export
cgr_sim <- function(g,
                    lvls        = NULL, # NULL if pure_cont
                    hom         = FALSE,
                    ns          = 10000,
                    nc          = 1,                    
                    alpha_rate = 10,
                    beta_rate  = 10,
                    sigma_rate = 10,                    
                    cell_rate  = 5,
                    beta_hom   = 10L,
                    sigma_hom  = 10L
                    ) {

  stopifnot(
    alpha_rate > 0,
    beta_rate  > 0,
    sigma_rate > 0,
    cell_rate  > 0
    )
  
  vd <- intersect(names(lvls), names(g))
  vc <- setdiff(names(g), vd)

  # Vertex info
  vs <- structure(list(nodes = c(vd, vc), disc = vd, cont = vc), class = c("verts", "list"))
  tv <- type_of_verts(vs) # if(!neq_empt_chr(vd)) "pure_cont" else if(!neq_empt_chr(vc)) "pure_cont" else "mixed"

  if (tv == "mixed") {
    if(!ess::is_decomposable(g, vd)) stop("g is not decomposable!")
  } else {
    if(!ess::is_decomposable(g)) stop("g is not decomposable!")
  }
  
  # Discrete part
  Ad_mat <- if (tv != "pure_cont") {
    gd <- if (neq_empt_chr(vc)) ess::subgraph(vc, g) else g
    sim_pure_disc_data(gd, lvls, ns, cell_rate)
  } else {NULL}
  
  # Parents 
  pa <- marginals_and_parents(g, tv, vs, Ad_mat)$pa

  # Pure cont case
  if (tv == "pure_cont") return(sim_pure_cont_data(pa, ns, coefs_pure_cont(pa, alpha_rate, beta_rate, sigma_rate)))

  # Pure disc case
  Ad <- as.data.frame(Ad_mat, stringsAsFactors = FALSE)
  if (tv == "pure_disc") return(Ad)

  # Mixed case
  coefs <- coefs_mixed(attr(pa, "nd"), hom, alpha_rate, beta_rate, sigma_rate, beta_hom, sigma_hom)
  Ac    <- sim_mixed_data(Ad, pa, ns, coefs, hom)
  sims  <- cbind(Ac, Ad)
  return(structure(sims, "pa" = pa, "coefs" = coefs))
}



## library(dplyr)
## library(ggplot2)
## library(stringi) # just for the nice %s+% infix function
## library(glue)

## g = list(
##   A = c("B", "X", "Y"),
##   B = c("A", "Y"),
##   X = c("A", "Y"),
##   Y = c("A", "X", "B")
## )

## d <- data.frame(A = "", B = "", X = 1, Y = 1)

## gp <- gengraph(d, "gen", g)

## plot(gp)

## lvls <- list(
##   A = c("0", "1"),# , "2", "3"),
##   B = c("0", "1") #, "2", "3")
## )

 
## set.seed(7)
## d1 <- cgr_sim(g, lvls, ns = 150,
##   alpha_rate = 10,
##   beta_rate  = 10,
##   sigma_rate = 10,
##   cell_rate  = 5) %>%
##   as_tibble()


## ggplot(d1, aes(x = Y, y = X)) +
##   geom_point(size = 0.2) +
##   geom_smooth(color = "red", method = "lm", se = FALSE) +
##   facet_wrap(~ A, scales = "free", nrow = 1, labeller = label_both) + 
##   theme_bw() +
##   theme(strip.background = element_rect(fill = "white"))
