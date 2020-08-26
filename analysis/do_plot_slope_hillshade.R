library(dplyr)
library(ggplot2)

d <- readRDS("../data/prc/covtype.Rds")

make_pair_data <- function(data, i, j) {
  .d <- data  %>%
    as_tibble() %>%
    filter(X55 %in% c({{i}}, {{j}})) %>%
    mutate(X55 = case_when(
      X55 == {{i}} ~ "black",
      X55 == {{j}} ~ "grey"
    ))
  colnames(.d)[1:3] <- c("Elevation", "Aspect", "Slope")
  .d
}

d12 <- make_pair_data(d, "1", "4")

ggplot(d12, aes(x = Elevation, y = Slope, color = X55)) +
  geom_point(size = 2) +
  scale_colour_identity() +
  theme_bw()
