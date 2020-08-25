verts <- function(d) {
  # d : data
  allowed_disc_class <- c("factor", "character") # remove factor when gRapHD is obsolete
  nodes <- colnames(d)
  node_class <- .map_chr(d, class)
  stopifnot(length(unique(node_class)) < 3L) # Either factor or character but not both!
  disc_idx <- which(node_class %in% allowed_disc_class)
  disc <- nodes[disc_idx]
  cont <- nodes %[int% -disc_idx
  lvls <- .map_int(seq_along(node_class), function(x) {
    if (node_class[x] %in% allowed_disc_class) {
      return(length(unique(d[, x, drop = TRUE])))
    }
    else {
      return(0L)
    }
  })
  structure(list(
    nodes = nodes,
    disc  = disc, # character(0) if empty
    cont  = cont, # character(0) if empty
    lvls  = structure(lvls, names = nodes)
  ),
  class = c("verts", "list")
  )
}

nodes <- function(v) UseMethod("nodes")
disc  <- function(v) UseMethod("disc")
cont  <- function(v) UseMethod("cont")
lvls  <- function(v) UseMethod("lvls")
nodes.verts <- function(v) v$nodes
disc.verts  <- function(v) v$disc
cont.verts  <- function(v) v$cont
lvls.verts  <- function(v) v$lvls

type_of_verts <- function(vs) {
  # vs: verts object
  vd <- disc(vs)
  vc <- cont(vs)
  ## if(!neq_empt_chr(vd)) return("pure_cont") else if(!neq_empt_chr(vc)) return("pure_cont") else return("mixed")
  lv <- lvls(vs)
  if (min(lv) == 0L && max(lv) > 0L) return("mixed")
  if (max(lv) == 0L) return("pure_cont")
  if (min(lv)  > 0L) return("pure_disc")
  stop("Something went wrong; data has neither discrete or continuous vertices!", call. = FALSE)
}
