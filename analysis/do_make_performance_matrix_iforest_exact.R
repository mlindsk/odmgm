library(IsolationForest)
library(dplyr)
library(glue)
library(stringi)

set.seed(300718)

classes <- readRDS("objects/classes.Rds")
nclass <- length(classes)
pm <- matrix(0L, nclass, nclass, dimnames = list(classes, classes))

for (row in classes) {

  d_row <- glue("objects/data_{row}.Rds") %>% readRDS()
  
  for (col in classes) {

    msg <- glue(" (row, col) = ({row}, {col})")
    cat(msg, "\n")
    
    d_col <- glue("objects/data_{col}.Rds") %>% readRDS()

    row_if    <- IsolationTrees(d_row, ntree = 100)
    row_score <- AnomalyScore(d_row, row_if)$outF
    col_score <- AnomalyScore(d_col, row_if)$outF
    crit_val  <- quantile(row_score, c(0, 0.95))[2]
    
    pm[row, col] <- mean(col_score >= crit_val)

    saveRDS(pm, "objects/pm_iforest_exact.Rds")
    
  }
}
