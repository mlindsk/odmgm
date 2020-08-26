library(molicMixed)
library(dplyr)
library(glue)

set.seed(300718)

classes <- readRDS("objects/classes.Rds")
nclass  <- length(classes)
ns <- 10000
nc <- 40

d <- readRDS("objects/data_full.Rds") %>%
  select(-X55)

## sorted_cols <- paste("X",
##   colnames(d) %>%
##     stringr::str_remove("X") %>%
##     as.integer() %>%
##     sort(), sep = "")


## d <- d %>% select(all_of(sorted_cols))

n_falses <- rep(FALSE, nrow(d))

d_outlier_in_which_class <- tibble(
  c1 = n_falses,
  c2 = n_falses,
  c3 = n_falses,
  c4 = n_falses,
  c5 = n_falses,
  c6 = n_falses,
  c7 = n_falses
)

count_try_errors <- 0L

for (i in 1:nrow(d)) {
  zi <- d[i, ]

  for (class in classes) {

    msg <- glue(" (i, class) = ({i}, {class}) - errors: {count_try_errors}")
    cat(msg, "\n")

    class_dat   <- paste(paste("objects/data_", class, sep = ""), ".Rds", sep = "")
    class_graph <- paste(paste("objects/graph_", class, sep = ""), ".Rds", sep = "")
    d_class <- readRDS(class_dat)
    g_class <- readRDS(class_graph)
    class_var <- paste("c", class, sep = "")

    m <- try(fit_outlier(d_class, zi, g_class, nsim = ns, ncores = nc, validate = FALSE))

    if (inherits(m, "try-error")) {
      count_try_errors <- count_try_errors + 1L
      d_outlier_in_which_class[i, class_var] <- FALSE # SHOULD BE NA INSTEAD! ?
    } else {
      d_outlier_in_which_class[i, class_var] <- m$pval <= 0.05
    }

    saveRDS(d_outlier_in_which_class, "objects/d_outlier_in_which_class.Rds")
    
  }
  
}
