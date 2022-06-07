#' Confounder blanket learner
#'
#' This function performs the confounder blanket learner (CBL) algorithm for
#' causal discovery.
#'
#' @param x Matrix or data frame of foreground variables. Currently only 
#'   supports continuous or binary features.
#' @param z Matrix or data frame of background variables. Currently only 
#'   supports continuous or binary features.
#' @param s Feature selection method. Includes native support for sparse linear
#'   regression (\code{s = "lasso"}) and gradient boosting (\code{s = "boost"}).
#'   Alternatively, a user-supplied function mapping features \code{x} and 
#'   outcome \code{y} to a bit vector indicating which features are selected. 
#'   See Examples.
#' @param B Number of complementary pairs to draw for stability selection. 
#'   Following Shah & Samworth (2013), we recommend leaving this fixed at 50.
#' @param gamma Omission threshold. If either of two foreground variables is 
#'   omitted from the model for the other with frequency \code{gamma} or higher,
#'   we infer that they are causally disconnected.
#' @param maxiter Maximum number of iterations to loop through if convergence
#'   is elusive. 
#' @param params Optional list to pass to \code{lgb.train} if \code{s = "boost"}.
#'   See \code{lightgbm::\link[lightgbm]{lgb.train}}.
#' @param parallel Compute stability selection subroutine in parallel? Must 
#'   register backend beforehand, e.g. via \code{doMC}. See example below.
#' @param ... Extra parameters to be passed to the feature selection subroutine.
#'   
#' @details 
#' The CBL algorithm (Watson & Silva, 2022) learns a partial order over 
#' foreground variables \code{x} via relations of minimal conditional 
#' (in)dependence with respect to a set of background variables \code{z}. The
#' method is sound and complete with respect to a so-called "lazy oracle", who 
#' only answers independence queries about variable pairs conditioned on the 
#' intersection of their respective non-descendants.
#' 
#' For computational tractability, CBL performs conditional independence tests 
#' via supervised learning with feature selection. The current implementation 
#' includes support for sparse linear models (\code{s = "lasso"}) and gradient 
#' boosting machines (\code{s = "boost"}). For statistical inference, CBL uses 
#' complementary pairs stability selection (Shah & Samworth, 2013), which bounds
#' the probability of errors of commission. 
#' 
#' @return 
#' A square, lower triangular ancestrality matrix. Call this matrix \code{m}. 
#' If CBL infers that $X_i \prec X_j$, then \code{m[j, i] = 1}. If CBL infers 
#' that $X_i \preceq X_j$, then \code{m[j, i] = 0.5}. If CBL infers that 
#' $X_i \sim X_j$, then \code{m[j, i] = 0}. Otherwise, \code{m[j, i] = NA}.
#' 
#' @references   
#' Watson, D.S. & Silva, R. (2022). Causal discovery under a confounder blanket.
#' \emph{arXiv} preprint, 2205.05715. 
#' 
#' Shah, R. & Samworth, R. (2013). Variable selection with error control: 
#' Another look at stability selection. \emph{J. R. Statist. Soc. B}, 
#' \emph{75}(1):55–80, 2013.
#' 
#' @examples 
#' # Load data
#' data(cbl_sim)
#' 
#' # Run CBL
#' m <- cbl(x, z)
#' 
#' # Run CBL in parallel
#' require(doMC)
#' registerDoMC(2)
#' m <- cbl(x, z, parallel = TRUE)
#' 
#' # With user-supplied feature selection subroutine
#' s_new <- function(x, y) {
#'   # Fit model, extract coefficients
#'   df <- data.frame(x, y)
#'   f_full <- lm(y ~ 0 + ., data = df)
#'   f_reduced <- step(f_full, trace = 0)
#'   keep <- names(coef(f_reduced))
#'   # Return bit vector
#'   out <- ifelse(colnames(x) %in% keep, 1, 0)
#'   return(out)
#' }
#' m <- cbl(x, z, s = s_new)
#' 
#' @export 
#' @import data.table
#' @import foreach

cbl <- function(
  x, 
  z, 
  s = 'lasso', 
  B = 50, 
  gamma = 0.5, 
  maxiter = NULL, 
  params = NULL,
  parallel = FALSE,
  ...
) {
  # Prelimz
  if (gamma < 0 | gamma > 1) {
    stop('gamma must be on [0,1].')
  }
  if (is.function(s)) {
    test_s <- s(z, x[, 1], ...)
    if (length(test_s) != ncol(z) | !all(test_s %in% c(0, 1))) {
      stop('s must provide a bit vector of length equal to the number of input features.')
    }
  }
  if (nrow(x) != nrow(z)) {
    stop('x and z must be the same sample size.')
  }
  if (is.data.frame(x)) {
    warning('Converting x to matrix format.')
    x <- model.matrix(~., x)[, -1]
  }
  if (is.data.frame(z)) {
    warning('Converting z to matrix format.')
    z <- model.matrix(~., z)[, -1]
  }
  if (B < 50) {
    warning('Results may be unstable with B < 50.')
  } else if (B > 50) {
    warning('The r-concavity assumption may not hold for B > 50.')
  }
  if (is.null(maxiter)) {
    maxiter <- Inf 
  }
  n <- nrow(x)
  d_x <- ncol(x)
  d_z <- ncol(z)
  
  # Initialize
  adj_list <- list(
    matrix(NA_real_, nrow = d_x, ncol = d_x, 
           dimnames = list(colnames(x), colnames(x)))
  )
  converged <- FALSE
  iter <- 0
  ### BIG LOOP ###
  while(converged == FALSE & iter <= maxiter) {
    # Extract relevant adjacency matrices
    if (iter == 0) {
      adj0 <- adj1 <- adj_list[[1]]
    } else {
      adj0 <- adj_list[[iter]]
      adj1 <- adj_list[[iter + 1]]
    }
    # Subsampling loop
    sub_loop <- function(b, i, j, a1) {
      z_t <- cbind(z, x[, a1])
      d_zt <- ncol(z_t)
      # Take complementary subsets
      a_set <- sample(n, round(0.5 * n))
      b_set <- seq_len(n)[-a_set]
      # Fit reduced models
      s0 <- sapply(c(i, j), function(k) {
        c(l0(z_t[a_set, ], x[a_set, k], s, params, ...), 
          l0(z_t[b_set, ], x[b_set, k], s, params, ...))
      })
      # Fit expanded models
      s1 <- sapply(c(i, j), function(k) {
        not_k <- setdiff(c(i, j), k)
        c(l0(cbind(z_t[a_set, ], x[a_set, not_k]), x[a_set, k], s, params, ...),
          l0(cbind(z_t[b_set, ], x[b_set, not_k]), x[b_set, k], s, params, ...))
      })
      # Record disconnections and (de)activations
      dis_a <- any(s1[d_zt + 1, ] == 0)
      dis_b <- any(s1[2 * (d_zt + 1), ] == 0)
      dis <- rep(c(dis_a, dis_b), each = d_zt)
      d_ji <- s0[, 1] == 1 & s1[seq_len(d_zt), 1] == 0
      a_ji <- s0[, 2] == 0 & s1[seq_len(d_zt), 2] == 1
      d_ij <- s0[, 2] == 1 & s1[seq_len(d_zt), 2] == 0
      a_ij <- s0[, 1] == 0 & s1[seq_len(d_zt), 1] == 1
      # Export
      out <- data.table(b = rep(c(2 * b - 1, 2 * b), each = d_zt), i, j,
                        z = rep(colnames(z_t), times = 2),
                        dis, d_ji, a_ji, d_ij, a_ij)
      return(out)
    }
    # Pairwise test loop
    for (i in 2:d_x) {
      for (j in 1:(i - 1)) {
        # Only continue if relationship is unknown
        if (is.na(adj1[i, j]) & is.na(adj1[j, i])) { 
          preceq_i <- which(adj0[i, ] > 0)
          preceq_j <- which(adj0[j, ] > 0)
          a0 <- intersect(preceq_i, preceq_j) 
          preceq_i <- which(adj1[i, ] > 0)
          preceq_j <- which(adj1[j, ] > 0)
          a1 <- intersect(preceq_i, preceq_j) 
          # Only continue if the set of non-descendants has increased since last 
          # iteration (i.e., have we learned anything new?)
          if (iter == 0 | length(a1) > length(a0)) {
            if (parallel == TRUE) {
              df <- foreach(bb = seq_len(B), .combine = rbind) %dopar%
                sub_loop(bb, i, j, a1)
            } else {
              df <- foreach(bb = seq_len(B), .combine = rbind) %do%
                sub_loop(bb, i, j, a1)
            }
            # Compute rates
            df[, disr := sum(dis) / .N]
            if (df$disr[1] > gamma) { 
              adj1[i, j] <- adj1[j, i] <- 0
            } else {
              df[, drji := sum(d_ji) / .N, by = z]
              df[, arji := sum(a_ji) / .N, by = z]
              df[, drij := sum(d_ij) / .N, by = z]
              df[, arij := sum(a_ij) / .N, by = z]
              df <- unique(df[, .(i, j, z, disr, drji, arji, drij, arij)])
              # Consistent lower bound
              lb <- lb_fn(df, B)
              # Stable upper bound
              out <- foreach(oo = c('ji', 'ij'), .combine = rbind) %:%
                foreach(rr = c('R1', 'R2'), .combine = rbind) %do%
                ss_fn(df, lb, oo, rr, B)
              # Update adjacency matrix
              if (sum(out$decision) == 1) {
                if (out[decision == 1, order == 'ji' & rule == 'R1']) {
                  adj1[i, j] <- 1
                  adj1[j, i] <- 0
                } else if (out[decision == 1, order == 'ji' & rule == 'R2']) {
                  adj1[i, j] <- 0.5
                  adj1[j, i] <- 0
                } else if (out[decision == 1, order == 'ij' & rule == 'R1']) {
                  adj1[j, i] <- 1
                  adj1[i, j] <- 0
                } else if (out[decision == 1, order == 'ij' & rule == 'R2']) {
                  adj1[j, i] <- 0.5
                  adj1[i, j] <- 0
                }
              } else if (sum(out$decision == 2)) {
                if (out[order == 'ji', sum(decision) == 2]) {
                  adj1[i, j] <- 0.5
                } else if (out[order == 'ij', sum(decision) == 2]) {
                  adj1[j, i] <- 0.5
                }
              }
            }
          }
        } 
      }
    }
    # Store that iteration's adjacency matrix
    iter <- iter + 1
    adj_list <- append(adj_list, list(adj1))
    # Check for convergence
    if (identical(adj0, adj1)) {
      converged <- TRUE
    }
  }
  # Export final adjacency matrix
  adj_mat <- adj_list[[length(adj_list)]]
  return(adj_mat)
}


