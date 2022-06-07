#' CPSS utility functions
#' 
#' Compute the tail probability of an r-concave random variable. Code taken
#' verbatim from Rajen Shah's personal website:
#' \url{http://www.statslab.cam.ac.uk/~rds37/papers/r_concave_tail.R}.
#' 
#' @param eta Upper bound on the expectation of the r-concave random variable.
#' @param B Number of complementary pairs for stability selection. 
#' @param r Of r-concavity fame.
#'
#' @importFrom stats optimize uniroot

r.TailProbs <- function(eta, B, r) {
  # TailProbs returns a vector with the tail probability for each \tau = ceil{B*2\eta}/B + 1/B,...,1
  # We return 1 for all \tau = 0, 1/B, ... , ceil{B*2\eta}/B
  MAXa <- 100000
  MINa <- 0.00001
  s <- -1/r
  etaB <- eta * B
  k_start <- (ceiling(2 * etaB) + 1)
  output <- rep(1, B)
  if(k_start > B) return(output)
  
  a_vec <- rep(MAXa,B)
  
  Find.a <- function(prev_a) uniroot(Calc.a, lower = MINa, upper = prev_a, tol = .Machine$double.eps^0.75)$root
  Calc.a <- function(a) {
    denom <- sum((a + 0:k)^(-s))
    num <- sum((0:k) * (a + 0:k)^(-s))
    num / denom - etaB
  }
  
  for(k in k_start:B) a_vec[k] <- Find.a(a_vec[k-1])
  
  OptimInt <- function(a) {
    num <- (k + 1 - etaB) * sum((a + 0:(t-1))^(-s))
    denom <- sum((k + 1 - (0:k)) * (a + 0:k)^(-s))
    1 - num / denom
  }
  # NB this function makes use of several gloabl variables
  
  prev_k <- k_start
  for(t in k_start:B) {
    cur_optim <- rep(0, B)
    cur_optim[B] <- OptimInt(a_vec[B])
    if (prev_k <= (B-1)) {
      for (k in prev_k:(B-1))
        cur_optim[k] <- optimize(f=OptimInt,lower = a_vec[k+1], upper = a_vec[k], maximum  = TRUE)$objective
    }
    output[t] <- max(cur_optim)
    prev_k <- which.max(cur_optim)
  }
  return(output)
}

#' CPSS upper bound
#' 
#' Compute the min-D factor of Shah & Samworth's Eq. 8 (2013). Code taken
#' verbatim from Rajen Shah's personal website:
#' \url{http://www.statslab.cam.ac.uk/~rds37/papers/r_concave_tail.R}.
#' 
#' @param theta Low rate threshold.
#' @param B Number of complementary pairs for stability selection. 
#' @param r Of r-concavity fame.

minD <- function(theta, B, r = c(-1/2, -1/4)) {
  pmin(c(rep(1, B), r.TailProbs(theta^2, B, r[1])), r.TailProbs(theta, 2*B, r[2]))
}

