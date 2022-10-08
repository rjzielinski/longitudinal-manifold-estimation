norm_euclidean <- function(x) {
  # Norm function in a Euclidean space of any dimension
  # return(norm(matrix(x, ncol = 1), type = "F"))
  return(sqrt(sum(as.vector(x) ^ 2)))
}