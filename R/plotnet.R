#' Plot the wavelet correlation network.
#' @param mat a association strength matrix.
#' @return
#' @export
#' network graph.
#' @examples
#' plotnet(mat = cor(Fauna.data[,-12]))
plotnet <- function(mat){
  map <- function(x, range = c(0,1), from.range=NA) {
    if(any(is.na(from.range))) from.range <- range(x, na.rm=TRUE)
    if(!diff(from.range)) return(
      matrix(mean(range), ncol=ncol(x), nrow=nrow(x), dimnames = dimnames(x)) )
    x <- (x-from.range[1]);  x <- x/diff(from.range)
    if(diff(from.range) == 0) x <- 0
    if (range[1]>range[2]) x <- 1-x
    x <- x*(abs(diff(range))) + min(range)
    x[x<min(range) | x>max(range)] <- NA
    x
  }
  group <- collor[c(2,2,2,4,5,5,5,3,3,3,3)]
  comm.net <- igraph::graph.adjacency(mat, weighted = T, mode = "undirected")
  cols <- c("#0066FF", "#D55E00")
  igraph::E(comm.net)$color <- ifelse(igraph::E(comm.net)$weight < 0, cols[1], cols[2])
  igraph::E(comm.net)$width <- abs(igraph::E(comm.net)$weight) * 2
  E(comm.net)$label <- NA
  E(comm.net)$curved <- 0.3
  igraph::V(comm.net)$size <- map(degree(comm.net),c(5,10))*2
  igraph::V(comm.net)$color <- group
  igraph::V(comm.net)$frame.color <- group
  la <- igraph::layout.circle(comm.net)
  plot(comm.net, layout = la, vertex.label = E(comm.net)$label)
  x = la[,1]*1.37 ; y = la[,2]*1.37
  angle = ifelse(atan(-(la[,1]/la[,2]))*(180/pi) < 0,
                 90 + atan(- (la[,1]/la[,2]))*(180 / pi),
                 270 + atan(-la[,1]/la[,2])*(180 / pi))
  for (i in 1:length(x)) {
    text(x = x[i], y = y[i], labels = rownames(mat)[i], col = "black",
         cex = 1.4, srt = angle[i], xpd = T)
  }
}
