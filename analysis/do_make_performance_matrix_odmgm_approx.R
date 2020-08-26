library(molicMixed)
library(dplyr)
library(glue)

classes <- readRDS("objects/classes.Rds")
nclass <- length(classes)
pm <- matrix(0L, nclass, nclass, dimnames = list(classes, classes))
ns <- 10000
nc <- 3

for (row in classes) {

  row_dat   <- paste(paste("objects/data_", row, sep = ""), ".Rds", sep = "")
  row_graph <- paste(paste("objects/graph_", row, sep = ""), ".Rds", sep = "")
  d_row <- readRDS(row_dat)
  g_row <- readRDS(row_graph)

  row_model <- fit_outlier_model(d_row, g_row, nsim = ns, ncores = nc, validate = FALSE)
  
  for (col in classes) {

    msg <- glue(" (row, col) = ({row}, {col})")
    cat(msg, "\n")
    
    col_dat <- paste(paste("objects/data_", col, sep = ""), ".Rds", sep = "")
    d_col   <- readRDS(col_dat)

    is_outlier <- sapply(1:nrow(d_col), function(i) {
      z  <- d_col[i, ]
      dz <- deviance(row_model, z)
      pv <- pval(row_model, dz)
      pv <= 0.05
    })

    pm[row, col] <- is_outlier %>% mean()
    
    saveRDS(pm, "objects/pm_approx.Rds")
  }
}
