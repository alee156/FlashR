#' Adjacency Spectral Embedding
#'
#' The eigendecomposition of an adjacency matrix provides a way to embed a graph  as points in finite dimensional Euclidean space.  This embedding allows the full arsenal of statistical and machine learning methodology  for multivariate Euclidean data to be deployed for graph inference. Adjacency Spectral Embedding performs this embedding of matrix A into dim dimensions.
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
ase <- function(A, dims){
  if(nrow(A) >= 400){
    require(irlba)
    A.svd <- irlba(A, nu = dims, nv = dims)
    A.svd.values <- A.svd$d[1:dims]
    A.svd.vectors <- A.svd$v[,1:dims]
    if(dim == 1)
      A.coords <- sqrt(A.svd.values) * A.svd.vectors
    else
      A.coords <- A.svd.vectors %*% diag(sqrt(A.svd.values))
  } else{
    A.svd <- svd(A)
    if(dims == 1)
      A.coords <- A.svd$v[,1] * sqrt(A.svd$d[1])
    else
      A.coords <- A.svd$v[,1:dims] %*% diag(sqrt(A.svd$d[1:dims]))
  }

  return(A.coords)
}
