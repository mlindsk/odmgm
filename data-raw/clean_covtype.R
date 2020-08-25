# https://archive.ics.uci.edu/ml/datasets/Covertype
library(dplyr)
library(readr)
library(usethis)

d_ <- "../inst/extdata/covtype.data.gz" %>%
  read_delim(col_names = FALSE, delim = ",")

set.seed(1234)
covtype <- d_ %>%
  slice(sample(1:581012, 20000)) %>%
  mutate_at(colnames(d_)[11:ncol(d_)], as.character)

use_data(covtype, overwrite = TRUE)
