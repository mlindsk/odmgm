## ---------------------------------------------------------
##                NON-EXPORTED HELPERS
## ---------------------------------------------------------
.sim_Q_cont <- function(pa, M) {
  sum(.map_dbl(names(pa), function(p) {
    paj_c <- attr(pa, "pc")[[p]]
    dfj       <- M - length(paj_c) - 1
    if (!(dfj > 0)) return(0L) # stop("Degrees of freedom is <= 0 .sim_Q_cont!")
    -M * log(stats::rbeta(1, dfj/2, 1/2))
  }))
}

.sim_Q_mixed_inhom <- function(pa, z) {
  sum(.map_dbl(names(pa), function(p) {
    n_paj     <- attr(pa, "nd")[[p]]  # If paj_d = character(0) -> n_a = nrow(A) as needed.
    paj_d     <- attr(pa, "pd")[[p]]
    paj_c     <- attr(pa, "pc")[[p]]
    n_paj_i0  <- if (neq_empt_chr(paj_d)) n_paj[paste(z[paj_d], collapse = "")] else n_paj
    dfj       <- n_paj_i0 - length(paj_c) - 1L
    if (!(dfj > 0L)) return(0L) # stop("Degrees of freedom is <= 0 in .sim_Q_mixed_inhom!")
    -n_paj_i0 * log(stats::rbeta(1, dfj/2, 1/2))
  }))
}

.sim_Q_mixed_hom <- function(pa) {
  N <- sum(attr(pa, "nd")[[1L]]) # All marginal totals equals the grand total
  sum(.map_dbl(names(pa), function(p) {
    n_paj     <- attr(pa, "nd")[[p]]
    paj_d     <- attr(pa, "pd")[[p]]
    paj_c     <- attr(pa, "pc")[[p]]
    dfj       <- N - length(n_paj) - length(paj_c)
    if (!(dfj > 0)) return(0L) # stop("Degrees of freedom is <= 0 in .sim_Q_mixed_hom!")
    -N * log(stats::rbeta(1, dfj/2, 1/2))
  }))
}

.sim_internal_mixed <- function(pa, zs, hom = FALSE) {
  if (hom) {
    return(.map_dbl(zs, function(z) .sim_Q_mixed_hom(pa))) # indep. of z!
  } else {
    return(.map_dbl(zs, function(z) .sim_Q_mixed_inhom(pa, z)))
  }
}

.sim_internal_pure_cont <- function(pa, M, nsim) {
  .map_dbl(1:nsim, function(i) .sim_Q_cont(pa, M))
}

.sim_internal_pure_disc <- function(A,           # character matrix
                                    cms,
                                    sms,
                                    type,        # (pure_disc, pure_disc_raw, .dgm_sim)
                                    nsim         = 10000,
                                    ncores       = 1) {
  y       <- replicate(nsim, vector("character", ncol(A)), simplify = FALSE)
  M       <- nrow(A)
  Hx_M    <- Hx_(M)
  Delta   <- colnames(A)
  C1_vars <- attr(cms[[1]], "vars")
  C1_idx  <- match(C1_vars, Delta)
  p_nC1   <- cms[[1]] / M
  yC1_sim <- sample(names(p_nC1), nsim, replace = TRUE, prob = p_nC1)
  if (!(length(cms) - 1L)) { # The complete graph
    yC1_sim <- lapply(strsplit(yC1_sim, ""), function(z) {names(z) = C1_vars; z})
    if (type == "pure_disc") yC1_sim <- 2 * (.map_dbl(yC1_sim, TY, cms, sms) - Hx_M) # D(Y) 
    return(yC1_sim)
  }
  doParallel::registerDoParallel(ncores)
  combine_ <- switch(type, "pure_disc"  = 'c', "pure_disc_raw" = 'list', ".dgm_sim" = 'rbind')
  y <- foreach::`%dopar%`(foreach::foreach(z = 1:nsim,
    .combine      = combine_,
    .multicombine = TRUE,
    .maxcombine   = nsim,
    .inorder      = FALSE), {
      y_sim_z <- y[[z]]
      y_sim_z[C1_idx] <- .split_chars(yC1_sim[1])
      for (k in 2:length(cms)) {
        nCk     <- cms[[k]]
        Ck_vars <- attr(nCk, "vars")     # Clique names
        Ck_idx  <- match(Ck_vars, Delta) # Where is Ck in Delta
        nSk     <- sms[[k]]              # For Sk = Ø we have that nSk = M
        Sk_vars <- attr(nSk, "vars")     # Separator names
        if (is.null(Sk_vars)) {          # For empty separators
          p_nCk_minus_nSk <- nCk / nSk   # nSk = M !
          y_sim_z[Ck_idx] <- .split_chars(sample(names(p_nCk_minus_nSk), 1L, prob = p_nCk_minus_nSk))
        } else {
          Sk_idx              <- match(Sk_vars, Delta)
          Sk_idx_in_Ck        <- match(Sk_vars, Ck_vars)
          Ck_idx_minus_Sk_idx <- Ck_idx[-Sk_idx_in_Ck]
          ySk                 <- y_sim_z[Sk_idx]
          nSk_ySk             <- na_ya(nSk, paste0(ySk, collapse = ""))
          nCk_given_Sk        <- n_b(nCk, structure(Sk_idx_in_Ck, names = ySk) )
          p_nCk_given_Sk_ySk  <- nCk_given_Sk / nSk_ySk # Cant be Inf, since ySk MUST be present since we simulated it
          y_sim_z[Ck_idx_minus_Sk_idx] <- .split_chars(sample(names(p_nCk_given_Sk_ySk), 1L, prob =  p_nCk_given_Sk_ySk))
        }
      }
      sim_z <- structure(y_sim_z, names = Delta)
      if (type == "pure_disc") {
        sim_z <- 2 * (TY(sim_z, cms, sms) - Hx_M)   # D(y) = 2*(T(y) - H(M))
      }
      sim_z
    })
  doParallel::stopImplicitCluster()
  return(y)
}

marginals_and_parents <- function(adj, type, vs, Ad = NULL) {
  .rip      <- if (type == "mixed") rip_marked(adj, disc(vs)) else rip(adj)
  .rip_disc <- if (type == "mixed") rip(ess::subgraph(cont(vs), adj)) else if (type == "pure_disc") .rip else NULL
  cms <- if (type != "pure_cont") a_marginals(Ad, .rip_disc$C) else NULL
  sms <- if (type != "pure_cont") a_marginals(Ad, .rip_disc$S) else NULL
  pa  <- if (type != "pure_disc") {
    p  <- .rip$PA[which(names(.rip$PA) %in% cont(vs))]  # .rip$PA[cont(vs)]
    pd <- lapply(p, function(x) intersect(x, disc(vs))) # Parent disc
    pc <- lapply(p, function(x) intersect(x, cont(vs))) # Parent cont
    nd <- if (type == "mixed") lapply(p, function(x) n_a(Ad[, intersect(x, disc(vs)), drop = FALSE])) else NULL # Parent disc tables
    structure(p, "pd" = pd, "pc" = pc, "nd" = nd)
  } else {
    NULL
  }
  return(list(cms = cms, sms = sms, pa = pa))
}

tmp_msg <- function() {
  message("Temporary msg: \n
1) Note that A is no longer A matrix; it must be a data.frame \n
2) Find a way to test if dfj > 0 in sim_Q_mixed!")
}

## ---------------------------------------------------------
##                   EXPORTED HELPERS
## ---------------------------------------------------------

#' Simulate observations from a discrete decomposable graphical model
#'
#' This function simulates observations from a discrete DGM using
#' data
#' 
#' @param A Character Matrix (data)
#' @param adj Adjacency list of a decomposable graph
#' @param nsim Number of simulations
#' @param ncores Number of cores to use in parallelization
#' @return This function returns a matrix of dimension \code{nsim x ncol(A)}
#' where each row correspond to a simulated observation from a DGM represented
#' by \code{adj}.
#' @export
dgm_sim <- function(A, adj, nsim = 10000, ncores = 1) {
  if (!is_matrix_chr(A)) stop("A is not a character matrix.")
  RIP   <- rip(adj)
  cms   <- a_marginals(A, RIP$C)
  sms   <- a_marginals(A, RIP$S)
  out   <- .sim_internal_pure_disc(A, cms, sms, type = ".dgm_sim", nsim = nsim, ncores = ncores)
  row.names(out) <- NULL
  return(out)
}

