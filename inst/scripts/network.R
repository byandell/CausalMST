suppressPackageStartupMessages({
  library(RColorBrewer)
  library(readr)
  library(dplyr)
  library(igraph)
})

nodes <- readr::read_csv(file.path("inst/doc", "quatrad_nodes.csv"))
links <- readr::read_csv(file.path("inst/doc", "quatrad_links.csv"))

net <- igraph::graph_from_data_frame(d=links, vertices=nodes, directed=T)
net <- igraph::simplify(net, remove.multiple = F, remove.loops = T)

jpeg("inst/doc/models.jpg", width = 720)
par(mfcol=c(1,1), mar=c(0,2,0,2))
plot(net,
     vertex.shape = "none",
     vertex.label = paste(names(igraph::V(net)),
                          igraph::V(net)$model, sep = ":"),
     vertex.label.font = 2,
     vertex.label.color = "black",
     vertex.label.cex = 1.5,
     edge.color = "gray60",
     layout = cbind(c(3,3,4,6,1,1,5,5,2,0) / 6,
                    c(5,1,3,3,6,0,6,0,3,3) / 6))
dev.off()

net_class <- function(links, nodes, quad = "mediator") {
  links <- links[!is.na(links[[quad]]),]
  ids <- unique(c(links$from, links$to))
  nodes <- nodes[nodes$id %in% ids, ]
  nodes <- nodes[order(nodes[[quad]]),]

  net <- igraph::graph_from_data_frame(d=links, vertices=nodes, directed=T)
  net <- igraph::simplify(net, remove.multiple = F, remove.loops = T)

  plot(net,
       vertex.shape = "none",
       vertex.label = igraph::V(net)$model,
       vertex.label.font = 2,
       vertex.label.color = "black",
       vertex.label.cex = 1.5,
       edge.color = "gray90",
       vertex.label.cex = igraph::V(net)$weight,
       edge.label = igraph::E(net)$class,
       edge.width = igraph::E(net)$weight,
       edge.label.cex = 2,
       layout = matrix(c(0.5,0,1,0.5,0,0.5,0.5,1),4,2))
  mtext(quad)
}

jpeg("inst/doc/quads.jpg", width=720)
par(mfcol=c(2,3), mar=c(1,1,2,1))
for(i in names(links)[5:10])
  net_class(links, nodes, i)
dev.off()
