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
  if (!is_decomposable(adj)) stop("The graph corresponding to adj is not decomposable!")
  if (!setequal(colnames(A), names(adj))) stop("Variables in A and the names of adj do not conform!")
  RIP   <- rip(adj)
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

#' Subgraph
#'
#' Construct a subgraph with a given set of nodes removed
#'
#' @param x Character vector of nodes to be removed
#' @param adj Adjacency list (named) or a neighbor matrix with dimnames given as the nodes
#' @examples
#' adj1 <- list(a = c("b", "d"), b = c("a", "c", "d"), c = c("b", "d"), d = c("a", "c", "b"))
#' # Toy data so we can plot the graph
#' d <- data.frame(a = "", b = "", c ="", d = "")
#' adj <- gengraph(d, type = "gen", adj = adj1)
#' plot(adj)
#' subgraph(c("c", "b"), adj1)
#' subgraph(c("b", "d"), as_adj_mat(adj1))
#' @export
subgraph <- function(x, adj) {
  # x: vector of nodes to delete
  if (inherits(adj, "matrix")) {
    keepers <- setdiff(dimnames(adj)[[1]], x)
    adj <- adj[keepers, keepers]
    return(adj)
  }
  else if (inherits(adj, "list")) {
    adj <- adj[-match(x, names(adj))]
    adj <- lapply(adj, function(e) {
      rm_idx <- as.vector(stats::na.omit(match(x, e)))
      if (neq_empt_int(rm_idx)) return(e[-rm_idx])
      return(e)
    })
    return(adj)
  }
  else {
    stop("adj must either be a matrix of an adjacency list.", call. = FALSE)
  }
}

#' Finds the components of a graph
#'
#' @param adj Adjacency list
#' @return A list with the elements being the components of the graph
#' @export
components <- function(adj) {
  nodes <- names(adj)
  comps <- list()
  comps[[1]] <- dfs(adj, nodes[1])
  while (TRUE) {
    new_comp  <- setdiff(nodes, unlist(comps))
    if (identical(new_comp, character(0))) return(comps)
    comps <- c(comps, list(dfs(adj[new_comp], new_comp[1])))
  }
  return(comps)
}

#' A test for decomposability in undirected (possibly marked) graphs
#'
#' This function returns \code{TRUE} if the graph is decomposable and \code{FALSE} otherwise
#'
#' @param adj Adjacency list of an undirected graph
#' @param disc_vars Character vector of discrete variables in adj if the graph is marked 
#' @examples
#' # Non-marked graphs:
#' # 4-cycle:
#' adj1 <- list(a = c("b", "d"), b = c("a", "c"), c = c("b", "d"), d = c("a", "c"))
#' is_decomposable(adj1) # FALSE
#' # Two triangles:
#' adj2 <- list(a = c("b", "d"), b = c("a", "c", "d"), c = c("b", "d"), d = c("a", "c", "b"))
#' is_decomposable(adj2) # TRUE
#'
#' # Marked graphs:
#' # A line graph with a forbidden path: the definition of non-decomposability in marked graphs
#' adj3 <- list(a = "b", b = c("a", "c"), c = "b")
#' is_decomposable(adj, c("a", "c"))
#' @export
is_decomposable <- function(adj, disc_vars = NULL) { 
  if (is.null(disc_vars)) {
    m <- try(mcs(adj, check = TRUE), silent = TRUE)
    if (inherits(m, "list")) return(TRUE)
      else return(FALSE)    
  } else {
    m <- try(mcs_marked(adj, disc_vars, check = TRUE), silent = TRUE)
    if (inherits(m, "list")) return(TRUE)
      else return(FALSE)    
  }
}

#' Make a complete graph
#'
#' A helper function to make an adjacency list corresponding to a complete graph
#'
#' @param nodes A character vector containing the nodes to be used in the graph
#' @examples
#' @export
make_complete_graph <- function(nodes) {
  structure(lapply(seq_along(nodes), function(k) {
    nodes[-which(nodes == nodes[k])]
  }), names = nodes)
}

#' Make a null graph
#'
#' A helper function to make an adjacency list corresponding to a null graph (no edges)
#'
#' @param nodes A character vector containing the nodes to be used in the graph
#' @examples
#' @export
make_null_graph <- function(nodes) {
  structure(lapply(seq_along(nodes), function(x) {
    character(0)
  }), names = nodes)
}
