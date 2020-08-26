library(dplyr)
library(molic)
library(ess)

d <- readRDS("../data/prc/covtype.Rds") %>%
  as_tibble() %>%
  select_if(is.character)

g <- fit_graph(d %>% filter(X55 == "7") %>% select(-X55))
plot(g)
