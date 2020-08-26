library(igraph)
library(dplyr)
library(purrr)
library(broom)
library(molicMixed)
library(stringi)
library(ggplot2)
source("helpers.R")

relative_path_to_fig <- "../../paper/incl/fig/"


## ---------------------------------------------------------
##                  FULL GRAPH
## ---------------------------------------------------------
d <- "../data/prc/covtype.Rds" %>% readRDS() %>% as_tibble()
g <- "objects/full_graph_with_class_variable.Rds" %>% readRDS()
g <- make_igraph(d, g, "black")
V(g)$color[55] <- "grey"

pdf(file = relative_path_to_fig %s+% "full_graph.pdf")
plot(g, layout = layout_as_star(g, center = "X55"), vertex.size = 10, vertex.label = NA)
title(parse(text = "G"), cex.main = 3, line = -27)
dev.off()

## ---------------------------------------------------------
##                CLASS DISTRIBUTION
## ---------------------------------------------------------
d %>%
  group_by(X55) %>%
  summarise(n = n(), prop = n() / nrow(.))

## ---------------------------------------------------------
##        INVESTIGATING DISCRETE VAR CONFIGURATIONS
## ---------------------------------------------------------
discrete_configurations <- function(d, class) {
  d %>%
    filter(X55 == class) %>%
    select_if(is_character) %>%
    select(-X55) %>%
    apply(1, function(x) paste0(x, collapse = "")) %>%
    table() %>% names() %>% stringr::str_replace_all("0", "\\.")
}

d1_conf <- discrete_configurations(d, "1"); d1_conf
d2_conf <- discrete_configurations(d, "2"); d2_conf
d3_conf <- discrete_configurations(d, "3"); d3_conf
d4_conf <- discrete_configurations(d, "4"); d4_conf
d5_conf <- discrete_configurations(d, "5"); d5_conf
d6_conf <- discrete_configurations(d, "6"); d6_conf
d7_conf <- discrete_configurations(d, "7"); d7_conf

intersect(d1_conf, d2_conf)
intersect(d1_conf, d3_conf)
intersect(d1_conf, d4_conf)
intersect(d1_conf, d5_conf)
intersect(d1_conf, d6_conf)
intersect(d1_conf, d7_conf)

intersect(d7_conf, d1_conf)
intersect(d7_conf, d2_conf)
intersect(d7_conf, d3_conf)
intersect(d7_conf, d4_conf)
intersect(d7_conf, d5_conf)
intersect(d7_conf, d6_conf)

## ---------------------------------------------------------
##                   CLASS 1
## ---------------------------------------------------------
d1 <- "objects/data_1.Rds"  %>% readRDS()
g1 <- "objects/graph_1.Rds" %>% readRDS()
g1 <- make_igraph(d1, g1, "black")

.order <- setdiff(names(V(g)), setdiff(names(V(g)), names(V(g1))))

pdf(file = relative_path_to_fig %s+% "graph1.pdf")
plot(g1, layout = layout_in_circle(g1, order = .order), vertex.size = 10, vertex.label = NA)
title(parse(text = "G[1]"), cex.main = 3, line = -27)
dev.off()


## ---------------------------------------------------------
##                   CLASS 2
## ---------------------------------------------------------
d2 <- "objects/data_2.Rds"  %>% readRDS()
g2 <- "objects/graph_2.Rds" %>% readRDS()
g2 <- make_igraph(d2, g2)
plot(g2, layout = layout_in_circle(g2, order = .order), vertex.size = 10)


## ---------------------------------------------------------
##                   CLASS 3
## ---------------------------------------------------------
d3 <- paste(paste("objects/data_", 3, sep = ""), ".Rds", sep = "") %>% readRDS()
g3 <- paste(paste("objects/graph_", 3, sep = ""), ".Rds", sep = "") %>% readRDS()
g3 <- make_igraph(d3, g3)
plot(g3, layout = layout_in_circle(g3, order = .order), vertex.size = 10)


## ---------------------------------------------------------
##                   CLASS 4
## ---------------------------------------------------------
d4 <- "objects/data_4.Rds" %>% readRDS()
g4 <- "objects/graph_4.Rds" %>% readRDS()
g4 <- make_igraph(d4, g4)
plot(g4, layout = layout_in_circle(g4, order = .order), vertex.size = 10)


par4 <- molicMixed:::rip_marked(g4, molicMixed:::verts(d4)$disc)$PA
par4[molicMixed:::verts(d4)$cont]

plot_density_by_factors(d4, X3, X4, X17, X31)
plot_density_by_factors(d4, X1, NULL, X17, X34)
plot_density_by_factors(d4, X3, X4, X17, X31)
plot_density_by_factors(d4, X4, X5, X17, X31)
plot_density_by_factors(d4, X1, X10, X17, X31)


# Not good
construct_regression_models(d4, X8 ~ X3, X17, X34) %>%
  pull(fit) %>%
  map( ~ glance(.x))

# Ok
construct_regression_models(d4, X6 ~ X3 + X8, X17, X31) %>%
  pull(fit) %>%
  map( ~ glance(.x))

# Good
construct_regression_models(d4, X2 ~ X3 + X7 + X9, X17, X31) %>%
  pull(fit) %>%
  map( ~ glance(.x))

# Bad
construct_regression_models(d4, X10 ~ X1, X17, X31) %>%
  pull(fit) %>%
  map( ~ glance(.x))


## ---------------------------------------------------------
##                   CLASS 5
## ---------------------------------------------------------
d5 <- "objects/data_5.Rds" %>% readRDS()
g5 <- "objects/graph_5.Rds" %>% readRDS()
g5 <- make_igraph(d5, g5)
plot(g5, layout = layout_in_circle(g5, order = .order), vertex.size = 10)

par5 <- molicMixed:::rip_marked(g5, molicMixed:::verts(d5)$disc)$PA
par5[molicMixed:::verts(d5)$cont]

plot_density_by_factors(d5, X1, NULL, X13, X44)
plot_density_by_factors(d5, X4, NULL, X37, X44)


# This one is spot on
construct_regression_models(d5, X9 ~ X2 + X3 + X7 + X8, X13, X44) %>%
  pull(fit) %>%
  map( ~ summary(.x))


# Also quite good
construct_regression_models(d5, X2 ~ X3 + X8, X13, X44) %>%
  pull(fit) %>%
  map( ~ summary(.x))

# Bad
construct_regression_models(d5, X5 ~ X4, X37, X44) %>%
  pull(fit) %>%
  map( ~ summary(.x))


## ---------------------------------------------------------
##                   CLASS 6
## ---------------------------------------------------------
d6 <- "objects/data_6.Rds" %>% readRDS()
g6 <- "objects/graph_6.Rds" %>% readRDS()
g6 <- make_igraph(d6, g6)
plot(g6, layout = layout_in_circle(g6, order = .order), vertex.size = 10)

par6 <- molicMixed:::rip_marked(g6, molicMixed:::verts(d6)$disc)$PA
par6[molicMixed:::verts(d6)$cont]

plot_models(d6, X1, NULL, X14, X15, X20)
plot_models(d6, X5, NULL, X14, X24, X31)
plot_models(d6, X10, NULL, X13, X27)

# Ok.
construct_regression_models(d6, X4 ~ X5, X24, X31) %>%
  pull(fit) %>%
  map_dbl( ~ glance(.x) %>% pull(r.squared)) 

# Good!
construct_regression_models(d6, X8 ~ X3 + X7, X14, X20) %>%
  pull(fit) %>%
  map_dbl( ~ glance(.x) %>% pull(r.squared)) 

construct_regression_models(d6, X9 ~ X3 + X7 + X8, X14, X20) %>%
  pull(fit) %>%
  map( ~ glance(.x))

# ~ ok..!
construct_regression_models(d6, X7 ~ X3, X14, X20) %>%
  pull(fit) %>%
  map( ~ glance(.x))

# BAD!
construct_regression_models(d6, X3 ~ X1, X14, X20) %>%
  pull(fit) %>%
  map( ~ glance(.x))

# ok.
construct_regression_models(d6, X6 ~ X1, X14, X20) %>%
  pull(fit) %>%
  map( ~ glance(.x))

# Good!
construct_regression_models(d6, X2 ~ X3 + X8 + X9, X14) %>%
  pull(fit) %>%
  map( ~ glance(.x))


## ---------------------------------------------------------
##                   CLASS 7
## ---------------------------------------------------------
d7 <- "objects/data_7.Rds" %>% readRDS()
g7 <- "objects/graph_7.Rds" %>% readRDS()
g7 <- make_igraph(d7, g7, "black")
plot(g7, layout = layout_in_circle(g7, order = .order), vertex.size = 10, vertex.label = NA)

pdf(file = relative_path_to_fig %s+% "graph7.pdf")
plot(g7, layout = layout_in_circle(g7, order = .order), vertex.size = 10, vertex.label = NA)
title(parse(text = "G[7]"), cex.main = 3, line = -27)
dev.off()

## Subgraph

h7 <- induced_subgraph(g7,
  c("X1", "X2","X3","X5","X6","X7","X8","X9","X10","X11", "X12", "X52","X53","X54")
)


pdf(file = relative_path_to_fig %s+% "subgraph7.pdf")
plot(h7,
  vertex.size = 10,
  layout = layout_in_circle(h7),
  vertex.label.dist = 1.5,
  vertex.label.cex  = 2,
  vertex.label.color = "black",
  vertex.label = parse(text = paste0("v[", stringr::str_remove(V(h7)$name, "X"), "]")),
  label.family = "times"
)
title(parse(text = "H[7]"), cex.main = 2.5, line = -27)
dev.off()

# Plotting models

par7 <- molicMixed:::rip_marked(g7, molicMixed:::verts(d7)$disc)$PA
par7[molicMixed:::verts(d7)$cont]

plot_density_by_factors(d7, X1, NULL, X52, X54)
plot_density_by_factors(d7, X10, NULL, X11, X12)


# Bad
construct_regression_models(d7, X5 ~ X6 + X10, X11) %>%
  pull(fit) %>%
  map( ~ glance(.x))


# Ok
construct_regression_models(d7, X7 ~ X2 + X3, X52, X53) %>%
  pull(fit) %>%
  map( ~ glance(.x))

# Good!
construct_regression_models(d7, X8 ~ X2 + X3 + X7, X52, X53) %>%
  pull(fit) %>%
  map( ~ glance(.x))

# Perfect!
construct_regression_models(d7, X9 ~ X2 + X3 + X7 + X8, X52, X53) %>%
  pull(fit) %>%
  map( ~ glance(.x))


# Good!
construct_regression_models(d7, X6 ~ X10, X11, X12) %>%
  pull(fit) %>%
  map( ~ glance(.x))


# Bad!
construct_regression_models(d7, X3 ~ X2, X52, X53) %>%
  pull(fit) %>%
  map( ~ glance(.x))
