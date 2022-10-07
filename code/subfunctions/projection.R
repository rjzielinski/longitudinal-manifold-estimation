projection <- function(x, f, initial.guess) {
  DD <- function(t) {
    return(dist_euclidean(x, f(t)))
  }
  est <- nlm(DD, p = initial.guess)
  return(est$estimate)
}