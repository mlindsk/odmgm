library(dplyr)
library(odmgm)

d <- covtypes

## ---------------------------------------------------------
##                MAKE THE FULL GRAPH
## ---------------------------------------------------------
# g <- odmgm:::fit_mixed_graph(d) # takes forever...
# saveRDS(g, "objects/full_graph_with_class_variable.Rds")

v  <- molicMixed:::verts(d)
cont <- v$cont
disc <- v$disc
vcol <- structure(vector("character", length(c(cont, disc))), names = c(cont, disc))
vcol[1:length(cont)]  <- "lightsteelblue2"
vcol[(length(cont) + 1):ncol(d)]  <- "orange"
plot(gengraph(d, "gen", g), vcol, vertex.size = 2)

## ---------------------------------------------------------
##                   STORE NEW DATA       
## ---------------------------------------------------------
important_vars <- components(g)[[1]]
d <- d[, important_vars]
saveRDS(d, "objects/data_full.Rds")

## ---------------------------------------------------------
##                 MAKE THE CLASS GRAPHS
## ---------------------------------------------------------
classes <- d %>% select(X55) %>% pull() %>% unique() %>% sort()
saveRDS(classes, "objects/classes.Rds")

for (class in classes) {

  print(class)
  
  class_data <- d %>% filter(X55 == class) %>% select(-X55)
  file_name_data <- paste(paste("objects/data_", class, sep = ""), ".Rds", sep = "")
  saveRDS(class_data, file_name_data)

  class_graph     <- fit_mixed_graph(class_data)
  file_name_graph <- paste(paste("objects/graph_", class, sep = ""), ".Rds", sep = "")
  saveRDS(class_graph, file_name_graph)
  
}
