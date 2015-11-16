#' Compute Distance
#'
#' Compute the distance
# @seealso
#' @param A matrix, an adjacency matrix that used for eigendecomposition
#' @param dims dimensions for which the Euclidean space is defined
#' @return Adjacency Spectral Embedding
# @seealso
# @export
#' @examples
#' N <- 200
#' A <- matrix(data = runif(N*N), ncol = N, nrow = N)
#' dims <- 15
#' A_embed <- ase(A, dims)

compute_distance <- function(graphs, normx='F') {
  S <- dim(graphs)[3]
  dist <- matrix(rep(0, S*S), ncol=S)
  for (i in 1:dim(graphs)[3]) {
    for (j in i:dim(graphs)[3]) {
      dist[i,j] <- norm(graphs[,,i]-graphs[,,j], normx)
    }
  }
  dist <- dist + t(dist)
  return(dist)
}
