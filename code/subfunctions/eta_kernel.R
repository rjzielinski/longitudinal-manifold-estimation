eta_kernel <- function(t, lambda) {
  # Reproducing kernels associated with Sobolev space D^{-2}L^2(R^d)
  if (lambda %% 2 == 0) {
    if (norm_euclidean(t) == 0) {
      y <- 0
    } else {
      y <- (norm_euclidean(t) ^ lambda) * log(norm_euclidean(t))
    }
  } else {
    y <- norm_euclidean(t) ^ lambda
  }
  return(y)
  
  # norm_val <- norm_euclidean(t)
  # if (lambda %% 2 == 0) {
  #   if (norm_val == 0) {
  #     y <- 0
  #   } else {
  #     y <- (norm_val ^ lambda) * log(norm_val)
  #   }
  # } else {
  #   y <- norm_val ^ lambda
  # }
  # 
  # return(y)
}