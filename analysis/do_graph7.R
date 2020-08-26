## ---------------------------------------------------------
##                     LIBS
## ---------------------------------------------------------
library(igraph)
library(dplyr)
library(purrr)
library(broom)
library(molicMixed)
library(stringi)
library(ggplot2)
library(glue)
library(patchwork)
source("helpers.R")

relative_path_to_fig <- "../../paper/incl/fig/"

d7 <- "objects/data_7.Rds" %>% readRDS()
g7 <- "objects/graph_7.Rds" %>% readRDS()
g7 <- make_igraph(d7, g7, "black")


## ---------------------------------------------------------
##                 Gaussian densities
## ---------------------------------------------------------

d7_1 <- d7 %>%
  select(all_of(c("X1", "X52", "X54"))) %>%
  group_by_at(vars("X52", "X54")) %>%
  mutate(n = n()) %>%
  ungroup() %>%
  mutate(d12 = glue("list(v[52]=={X52}, v[54]=={X54},  n({X52}, {X54}) == {n})")) %>%
  select(-n)


p7_1 <- d7_1 %>%
  ggplot(aes(x = X1)) +
  geom_density() +
  facet_wrap(~ d12,
    scales = "free",
    strip.position = "top",
    labeller = label_parsed
  ) +
  theme_bw() +
  theme(strip.background =element_rect(fill = "white")) +
  theme(
    axis.text.y  = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  theme(strip.text = element_text(size = 12))


p7_1 <- p7_1 + theme(axis.title.y = element_blank(), axis.title = element_text(size = 12))
p7_1 <- p7_1 + xlab(parse(text = "v[1]")) 
p7_1

d7_10 <- d7 %>%
  select(all_of(c("X10", "X11", "X12"))) %>%
  group_by_at(vars("X11", "X12")) %>%
  mutate(n = n()) %>%
  ungroup() %>%
  mutate(d12 = glue("list(v[11]=={X11}, v[12]=={X12},  n({X11}, {X12}) == {n})")) %>%
select(-n)

p7_10 <- d7_10 %>%
  ggplot(aes(x = X10)) +
  geom_density() +
  facet_wrap(~ d12,
    scales = "free",
    strip.position = "top",
    labeller = label_parsed
  ) +
  theme_bw() +
  theme(strip.background =element_rect(fill = "white")) +
  theme(
    axis.text.y  = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  theme(strip.text = element_text(size = 12))

p7_10 <- p7_10 + theme(axis.title.y = element_blank(), axis.title = element_text(size = 12))
p7_10 <- p7_10 + xlab(parse(text = "v[10]"))
p7_10


p <- p7_1 / p7_10
p
ggsave(relative_path_to_fig %s+% "gaussian_given_discrete_parents.pdf", width = 7, height = 5)

# ?plotmath

## ---------------------------------------------------------
##             Multiple regressions
## ---------------------------------------------------------

# Bad
construct_regression_models(d7, X5 ~ X6 + X10, X11) %>%
  pull(fit) %>%
  map( ~ glance(.x) %>% pull(r.squared))

# Ok
construct_regression_models(d7, X7 ~ X2 + X3, X52, X53) %>%
  pull(fit) %>%
  map( ~ glance(.x) %>% pull(r.squared))

# Good!
construct_regression_models(d7, X8 ~ X2 + X3 + X7, X52, X53) %>%
  pull(fit) %>%
  map( ~ glance(.x) %>% pull(r.squared))

# Perfect!
construct_regression_models(d7, X9 ~ X2 + X3 + X7 + X8, X52, X53) %>%
  pull(fit) %>%
  map( ~ glance(.x) %>% pull(r.squared))
