# data_in_classes

library(UpSetR)
library(tidyverse)
library(odmgm)

d_cov <- covtype

d_class <- readRDS("objects/d_outlier_in_which_class.Rds") %>% 
  bind_cols(CT = paste(d_cov$X55)) %>% 
  select(CT, everything()) %>% 
  mutate_at(-1, ~as.integer(!.))

d_class %>% select(-CT) %>% colSums()

d_class %>% # select(-CT) %>% 
  distinct() %>% 
  arrange(CT)

tab_binary <- d_class %>% 
  unite(state, c1:c7, sep = "") %>% 
  count(state)

tab_sets <- d_class %>% 
  select(-CT) %>% 
  {(.)*col(.)} %>% 
  unite(state, c1:c7, sep = "") %>% 
  count(state) %>% 
  mutate(state = gsub("0", "", state))

tab_sets %>% mutate(set_size = nchar(state)) %>% 
  arrange(set_size, state)

tab_sets %>% pull(n) %>% sum()
  
upset(as.data.frame(d_class), nsets = 7, sets = paste0("c",7:1), keep.order = TRUE,
      order.by = c("freq", "degree"), group.by = "degree")

d_class_ <- d_class %>% 
  set_names(c("CT", paste("CT =", 1:7)))

source("counter_.R")
source("upset_.R")

upset_(
  as.data.frame(d_class_), nsets = 7, sets = paste("CT =",7:1), keep.order = TRUE,
  sets.x.label = "Samples not rejected in given CT",
  mainbar.y.label = "Sets of non-rejected CTs",
  set_size.show = TRUE, 
  point.size = 4, text.scale = 2.5, line.size = 1,
  set_size.scale_max = 19000,
  order.by = c("freq", "degree"), group.by = "degree"
)
dev.off()
