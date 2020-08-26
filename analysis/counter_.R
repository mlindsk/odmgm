Counter_ <- function (data, num_sets, start_col, name_of_sets, nintersections, 
          mbar_color, order_mat, aggregate, cut, empty_intersects, decrease){
  temp_data <- list()
  Freqs <- data.frame()
  end_col <- as.numeric(((start_col + num_sets) - 1))
  for (i in 1:num_sets) {
    temp_data[i] <- match(name_of_sets[i], colnames(data))
  }
  Freqs <- data.frame(plyr::count(data[, as.integer(temp_data)]))
  colnames(Freqs)[1:num_sets] <- name_of_sets
  if (is.null(empty_intersects) == F) {
    empty <- rep(list(c(0, 1)), times = num_sets)
    empty <- data.frame(expand.grid(empty))
    colnames(empty) <- name_of_sets
    empty$freq <- 0
    all <- rbind(Freqs, empty)
    Freqs <- data.frame(all[!duplicated(all[1:num_sets]), ], check.names = F)
  }
  ## Keep the "null" rows in the data
  ## Freqs <- Freqs[!(rowSums(Freqs[, 1:num_sets]) == 0), ]
  if (tolower(aggregate) == "degree") {
    for (i in 1:nrow(Freqs)) {
      Freqs$degree[i] <- rowSums(Freqs[i, 1:num_sets])
    }
    order_cols <- c()
    for (i in 1:length(order_mat)) {
      order_cols[i] <- match(order_mat[i], colnames(Freqs))
    }
    for (i in 1:length(order_cols)) {
      logic <- decrease[i]
      Freqs <- Freqs[order(Freqs[, order_cols[i]], decreasing = logic), ]
    }
  }
  else if (tolower(aggregate) == "sets") {
    Freqs <- Get_aggregates(Freqs, num_sets, order_mat, cut)
  }
  ## browser()
  delete_row <- (num_sets + 2)
  Freqs <- Freqs[, -delete_row]
  for (i in 1:nrow(Freqs)) {
    Freqs$x[i] <- i
    Freqs$color <- mbar_color
  }
  if (is.na(nintersections)) {
    nintersections = nrow(Freqs)
  }
  Freqs <- Freqs[1:nintersections, ]
  Freqs <- na.omit(Freqs)
  return(Freqs)
}