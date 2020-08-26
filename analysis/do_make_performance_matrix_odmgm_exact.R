library(molicMixed)
library(dplyr)
library(glue)

set.seed(300718)

classes <- readRDS("objects/classes.Rds")
nclass <- length(classes)
pm <- matrix(0L, nclass, nclass, dimnames = list(classes, classes))
ns <- 100
nc <- 1

for (row in classes) {

  row_dat   <- paste(paste("objects/data_", row, sep = ""), ".Rds", sep = "")
  row_graph <- paste(paste("objects/graph_", row, sep = ""), ".Rds", sep = "")
  d_row <- readRDS(row_dat)
  g_row <- readRDS(row_graph)
  
  for (col in classes) {

    col_dat <- paste(paste("objects/data_", col, sep = ""), ".Rds", sep = "")
    d_col   <- readRDS(col_dat)

    outlier_models <- sapply(1:nrow(d_col), function(i) {
      msg <- glue(" (row, col) = ({row}, {col})   -   {i}:{nrow(d_col)}")
      cat(msg, "\n")

      zi      <- d_col[i, ]
      d_row_i <- d_row
      if (row == col) d_row_i <- d_row_i[-i, ]

      m <- try(fit_outlier(d_row_i, zi, g_row, nsim = ns, ncores = nc, validate = FALSE))
      if (inherits(m, "try-error")) { # Some cells might not have enough information
        return(FALSE)
      } else {
        return(m$pval <= 0.05)
      }
    })

    pm[row, col] <- outlier_models %>% mean()

    saveRDS(pm, "objects/pm_exact.Rds")
    
  }
}
