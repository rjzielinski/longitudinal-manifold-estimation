ker <- function(x, mu, sigma) {
  # Smoothing kernel for density estimation
  yseq <- sapply((x - mu) / sigma, dnorm)
  return((sigma ^ {-length(x)}) * prod(yseq))
}