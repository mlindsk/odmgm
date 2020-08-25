gRapHD_empty_graph <- function(v, hom = FALSE) {
  # verts object
  methods::new("gRapHD",
    edges     = matrix(integer(0), , 2),
    p         = length(nodes(v)),
    numCat    = lvls(v),
    homog     = hom,
    vertNames = nodes(v)
  )
}

gRapHD_to_adj <- function(v, m) {
  # v: verts obj
  # m: gRapHD model
  A <- structure(gRapHD::adjMat(m), dimnames = list(nodes(v), nodes(v)))
  as_adj_lst(A)
}

#' Fit a marked decomposable graphical model
#' @description A generic method for structure learning in marked decomposable graphical models
#' @param d data.frame
#' @param tree Logical indicating whether or not the function should fit a tree or a general graph
#' @param hom Logical.
#' @param st Character. One of LR, AIC, or BIC
#' @return An adjacency list
#' @examples
#' \dontrun{
#' # TBA
#' }
#' @export
fit_mixed_graph <- function(d, tree = FALSE, hom = FALSE, st = "BIC") { # st in (LR, AIC, or BIC)
  if (!all(stats::complete.cases(d))) stop("data contains NA values.")
  v <- verts(d)
  if (inherits(d, "tbl_df")) d <- as.data.frame(d)
  d[] <- lapply(d, function(x) if(inherits(x, "character")) return(as.factor(x)) else return(x))

  m <- if (!tree) {
    gRapHD::stepw(
      model   = gRapHD_empty_graph(v, hom),
      dataset = as.data.frame(d),
      stat    = st,
      join    = TRUE
    )
  } else {
    gRapHD::minForest(d, homog = hom, stat = st)
  }
  return(gRapHD_to_adj(v, m))
}
