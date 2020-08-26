## ---------------------------------------------------------
##                EXPORTED HELPERS
## ---------------------------------------------------------
#' Probability mass function for a decomposable graphical model
#'
#' @param A A character matrix of data
#' @param adj Adjacency list or gengraph object of a decomposable graph. See \code{fit_graph}.
#' @param logp Logical; if TRUE, probabilities p are given as log(p).
#' @return A function - the probability mass function corresponding
#' to the decomposable graph \code{adj} using \code{A} to estimate
#' the probabilities.
#' @details It should be noted, that \code{pmf} is a closure function; i.e. its return
#' value is a function. Once \code{pmf} is called, one can query probabilities from the
#' returned function. If the probability of an observation is zero, the function will return
#' \code{NA} if \code{logp} is set to \code{TRUE} and zero otherwise.
#' @export
pmf <- function(A, adj, logp = FALSE) {
  if (!ess::is_decomposable(adj)) stop("The graph corresponding to adj is not decomposable!")
  if (!setequal(colnames(A), names(adj))) stop("Variables in A and the names of adj do not conform!")
  RIP   <- ess::rip(adj)
  cms   <- RIP$C
  sms   <- RIP$S
  ncms <- a_marginals(A, cms)
  nsms <- a_marginals(A, sms)
  .pmf <- function(y) {
    ny <- vapply(seq_along(ncms), FUN.VALUE = 1, FUN =  function(i) {
      nci    <- ncms[[i]]
      nsi    <- nsms[[i]]
      yci    <- y[match(attr(nci, "vars"), names(y))]
      ycinam <- paste0(yci, collapse = "")
      nciy   <- nci[ycinam] # NA if y is not seen on ci
      if (i == 1L) return(log(nciy))
      nsiy <- nsi[1]
      if (length(nsi) > 1) {
        ysi    <- y[match(attr(nsi, "vars"), names(y))]
        nsinam <- paste0(ysi, collapse = "")
        nsiy   <- nsi[nsinam] # NA if not seen on si
      } 
      return(log(nciy) - log(nsiy))
    })
    if (anyNA(ny)) return(ifelse(logp, NA, 0L)) # The observation was not seen on some marginals
    logprob <- sum(ny) - log(nrow(A))
    return(ifelse(logp, logprob , exp(logprob)))
  }
}
