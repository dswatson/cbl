#' Confounder blanket learner
#'
#' This function performs the confounder blanket learner (CBL) algorithm for
#' causal discovery.
#'
#' @param x Matrix or data frame of foreground variables. 
#' @param z Matrix or data frame of background variables. 
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
#'   register backend beforehand, e.g. via \code{doMC}. 
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
#' If CBL infers that \eqn{X_i \prec X_j}, then \code{m[j, i] = 1}. If CBL 
#' infers that \eqn{X_i \preceq X_j}, then \code{m[j, i] = 0.5}. If CBL infers 
#' that \eqn{X_i \sim X_j}, then \code{m[j, i] = 0}. Otherwise, 
#' \code{m[j, i] = NA}.
#' 
#' @references   
#' Watson, D.S. & Silva, R. (2022). Causal discovery under a confounder blanket.
#' To appear in \emph{Proceedings of the 38th Conference on Uncertainty in 
#' Artificial Intelligence}. \emph{arXiv} preprint, 2205.05715. 
#' 
#' Shah, R. & Samworth, R. (2013). Variable selection with error control: 
#' Another look at stability selection. \emph{J. R. Statist. Soc. B}, 
#' \emph{75}(1):55–80, 2013.
#' 
#' @examples 
#' # Load data
#' data(bipartite)
#' x <- bipartite$x
#' z <- bipartite$z
#' 
#' # Set seed
#' set.seed(42)
#' 
#' # Run CBL
#' cbl(x, z)
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
#' \donttest{
#' cbl(x, z, s = s_new)
#' }
#' 
#' @export 
#' @import data.table
#' @import foreach
#' @importFrom stats model.matrix

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
    #warning('Converting x to matrix format.')
    x <- model.matrix(~., x)[, -1]
  }
  if (is.data.frame(z)) {
    #warning('Converting z to matrix format.')
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
  # Nullify for visible bindings
  disr <- dis <- drji <- d_ji <- arji <- a_ji <- drij <- d_ij <- arij <- a_ij <- 
    bb <- oo <- rr <- rule <- decision <- NULL
  # Initialize
  adj_list <- list(
    matrix(NA_real_, nrow = d_x, ncol = d_x, 
           dimnames = list(colnames(x), colnames(x)))
  )
  adj0 <- adj1 <- adj_list[[1]]
  converged <- FALSE
  iter <- 0
  ### BIG LOOP ###
  while(converged == FALSE & iter <= maxiter) {
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
            z_t <- cbind(z, x[, a1])
            if (parallel == TRUE) {
              df <- foreach(bb = seq_len(B), .combine = rbind) %dopar%
                sub_loop(bb, i, j, x, z_t, s, params, ...)
            } else {
              df <- foreach(bb = seq_len(B), .combine = rbind) %do%
                sub_loop(bb, i, j, x, z_t, s, params, ...)
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
              eps <- epsilon_fn(df, B)
              # Stable upper bound
              out <- foreach(oo = c('ji', 'ij'), .combine = rbind) %:%
                foreach(rr = c('R1', 'R2'), .combine = rbind) %do%
                ss_fn(df, eps, oo, rr, B)
              # Update adjacency matrix
              if (sum(out$decision) > 0) {
                if (out[rule == 'R1' & order == 'ji', decision == 1]) {
                  adj1[i, j] <- 1
                  adj1[j, i] <- 0
                } else if (out[rule == 'R1' & order == 'ij', decision == 1]) {
                  adj1[j, i] <- 1
                  adj1[i, j] <- 0
                } else if (out[rule == 'R2' & order == 'ji', decision == 1]) {
                  adj1[i, j] <- 0.5
                  adj1[j, i] <- 0
                } else if (out[rule == 'R2' & order == 'ij', decision == 1]) {
                  adj1[j, i] <- 0.5
                  adj1[i, j] <- 0
                } 
              }
            }
          }
        } 
      }
    }
    # Iterate, check for convergence
    iter <- iter + 1
    if (identical(adj0, adj1)) {
      converged <- TRUE
    } else {
      adj_list <- append(adj_list, list(adj1))
      adj0 <- adj_list[[iter]]
      adj1 <- adj_list[[iter + 1]]
      adj_list[[iter]] <- 0
    }
  }
  # Export final adjacency matrix
  return(adj1)
}

if(getRversion() >= '2.15.1')  utils::globalVariables('.')
