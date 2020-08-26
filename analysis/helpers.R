## ---------------------------------------------------------
##                     HELPERS
## ---------------------------------------------------------

make_igraph <- function(data, graph, disc_col = "grey", cont_col = "white") {
  graph <- as_adj_mat(graph) %>%
    igraph::graph_from_adjacency_matrix("undirected")
  v           <- molicMixed:::verts(data)
  cont        <- v$cont
  graph_names <- names(V(graph))
  col <- rep(disc_col, ncol(data))
  col[which(graph_names %in% cont)] <- cont_col
  V(graph)$color <- col
  graph
}

## group_by_disc <- function(data, .c, .d1, .d2) {

##   ## d1_number <- stringr::str_remove(.d1, "X")
##   ## d2_number <- stringr::str_remove(.d2, "X")
  
##   data %>%
##     select(all_of(c(.c, .d1, .d2))) %>%
##     group_by_at(vars(.d1, .d2)) %>%
##     mutate(n = n()) %>%
##     ungroup() %>%
##     mutate(d12 = glue("n({d1_number},{d2_number}) = {n}")) %>%
##     select(-n)
## }


plot_density_by_factors <- function(data, .x, .y = NULL, .f1 = NULL, .f2 = NULL, .f3 = NULL) {

  is_null <- try(is.null(.y))
  
  if (!inherits(is_null, "try-error")) {
    data %>%
      ggplot(aes(x = {{.x}})) +
      facet_wrap(vars({{.f1}}, {{.f2}}, {{.f3}}), scale = "free_x") +
      geom_histogram() +
      theme_bw()
  } else {
    data %>%
      ggplot(aes(x = {{.x}}, y = {{.y}})) +
      geom_point() +
      geom_smooth(method = "lm") +
      facet_wrap(vars({{.f1}}, {{.f2}}, {{.f3}}), scale = "free_x") +
      theme_bw() 
  }
}

construct_regression_models <- function(.data, .formula, .f1 = NULL, .f2 = NULL) {
  .data %>%
    group_by({{.f1}}, {{.f2}}) %>%
    tidyr::nest() %>%
    mutate(
      fit = map(data, ~ lm(.formula, data = .x))
    ) 
}


tikz <- function (graph, layout) {

  if (class(layout) == "function") layout <- layout(graph)

  layout <- layout / max(abs(layout))

  cat("\\tikzset{\n")
  cat("\tnode/.style={circle, inner sep=1mm,minimum size=0.5cm, draw, thick, black, fill=red!20,text=black},\n")
  cat("\tnondirectional/.style={very thick,black},\n")
  cat("}\n")
  cat("\n")

  cat("\\begin{tikzpicture}[scale=5]\n")

  for (i in seq_along(V(graph))) {
    label <- V(graph)[i]$label
    # TODO: If node is cont. dont fill of node is discrete fill with black!
    msg <- glue(
      "\t\\node [node] ",
      "(v<<i>>) at ",
      "(<<layout[i,1]>>, <<layout[i,2]>>)",
      "\t{$v_{<<i>>}$};",
      .open  = "<<",
      .close = ">>"
    )
    cat(msg, "\n")
  }
  
  adj = get.adjacency(graph)

  for (i in 1:nrow(adj)) {
    for (j in 1:ncol(adj)) {
      if (adj[i, j] & j > i) {
        msg <- glue("\t\\path [nondirectional]",
          "(v<<i>>) edge (v<<j>>);",
          .open  = "<<",
          .close = ">>"
        )
        cat(msg, "\n")
      }
    }
  }
  cat("\\end{tikzpicture}\n")
}
