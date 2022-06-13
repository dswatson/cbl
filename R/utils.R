#' Feature selection subroutine
#' 
#' This function fits a potentially sparse supervised learning model and returns
#' a bit vector indicating which features were selected.
#' 
#' @param x Design matrix.
#' @param y Outcome vector.
#' @param s Regression method. Current options are \code{"lasso"} or
#'   \code{"boost"}.
#' @param params Optional list of parameters to use when \code{s = "boost"}.
#' @param ... Extra parameters to be passed to the feature selection subroutine.
#' 
#' @import glmnet
#' @import lightgbm

l0 <- function(x, y, s, params, ...) {
  if (is.function(s)) {
    out <- s(x, y, ...)
  } else {
    n <- nrow(x)
    trn <- sample(n, round(0.8 * n))
    tst <- seq_len(n)[-trn]
    if (s == 'lasso') {
      fit <- glmnet(x[trn, ], y[trn], intercept = FALSE, ...)
      y_hat <- predict.glmnet(fit, newx = x[tst, ], s = fit$lambda)
      eps <- y_hat - y[tst]
      mse <- colMeans(eps^2)
      betas <- coef.glmnet(fit, s = fit$lambda)[-1, which.min(mse)]
      out <- ifelse(betas == 0, 0, 1)
    } else if (s == 'boost') {
      d_trn <- lgb.Dataset(x[trn, ], label = y[trn])
      d_tst <- lgb.Dataset.create.valid(d_trn, x[tst, ], label = y[tst])
      fit <- lgb.train(params = params, data = d_trn, valids = list(tst = d_tst), 
                       verbose = 0, ...)
      vimp <- lgb.importance(fit)
      out <- as.numeric(colnames(x) %in% vimp$Feature)
    } 
  }
  return(out)
}


#' Complementary pairs subsampling loop
#' 
#' This function executes one loop of the model quartet for a given pair of 
#' foreground variables and records any disconnections and/or (de)activations.
#' 
#' @param b Subsample index.
#' @param i First foreground variable index.
#' @param j Second foreground variable index.
#' @param x Matrix of foreground variables.
#' @param z_t Intersection of iteration-t known non-descendants for foreground
#'   variables \code{i} and \code{j}.
#' @param s Regression method. Current options are \code{"lasso"} or
#'   \code{"boost"}.
#' @param params Optional list of parameters to use when \code{s = "boost"}.
#' @param ... Extra parameters to be passed to the feature selection subroutine.
#' 
#' @import data.table

sub_loop <- function(b, i, j, x, z_t, s, params, ...) {
  # Prelimz
  n <- nrow(x) 
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


#' Computer the consistency lower bound
#' 
#' @param df Table of (de)activation rates.
#' @param B Number of complementary pairs to draw for stability selection.
#' 
#' @import data.table
#' @import foreach

epsilon_fn <- function(df, B) {
  # Nullify 
  dji <- drji <- aji <- arji <- dij <- drij <- aij <- arij <- tau <- tt <-
    int_err <- ext_err <- NULL
  # Loop through thresholds in search of inconsistencies
  err_check <- function(tau) {
    # Inferences at this value of tau
    df[, dji := ifelse(drji >= tau, 1, 0)]
    df[, aji := ifelse(arji >= tau, 1, 0)]
    df[, dij := ifelse(drij >= tau, 1, 0)]
    df[, aij := ifelse(arij >= tau, 1, 0)]
    # Internal consistency (for a single Z)
    df[, int_err := ifelse(dji + aji + dij + aij > 1, 1, 0)]
    int_err <- ifelse(sum(df$int_err) > 0, 1, 0)
    # External consistency (across multiple Z's)
    if (df[, sum(dji) > 0] & df[, sum(dij + aij) > 0]) {
      ext_err <- 1
    } else if (df[, sum(dij) > 0] & df[, sum(dji + aji) > 0]) {
      ext_err <- 1
    } else {
      ext_err <- 0
    }
    # Export
    out <- data.table('tau' = tau, 'int_err' = int_err, 'ext_err' = ext_err)
    return(out)
  }
  err_df <- foreach(tt = seq_len(2 * B) / (2 * B), .combine = rbind) %do% 
    err_check(tt)
  # Compute minimal thresholds, export
  epsilon <- err_df[int_err == 0 & ext_err == 0, min(tau)]
  return(epsilon)
}


#' Infer causal direction using stability selection
#' 
#' @param df Table of (de)activation rates.
#' @param epsilon Consistency lower bound, as computed by \code{epsilon_fn}.
#' @param order Causal order of interest, either \code{"ij"} or \code{"ji"}.
#' @param rule Inference rule, either \code{"R1"} or \code{"R2"}.
#' @param B Number of complementary pairs to draw for stability selection.
#' 
#' @import data.table

ss_fn <- function(df, epsilon, order, rule, B) {
  # Nullify
  drji <- arji <- drij <- arij <- detected <- surplus <- err_bound <- NULL
  # Find the right rate
  if (order == 'ji' & rule == 'R1') {
    r <- df[, drji]
  } else if (order == 'ji' & rule == 'R2') {
    r <- df[, arji]
  } else if (order == 'ij' & rule == 'R1') {
    r <- df[, drij]
  } else if (order == 'ij' & rule == 'R2') {
    r <- df[, arij]
  } 
  # Stability selection parameters
  theta <- mean(r)
  ub <- minD(theta, B) * sum(r <= theta)
  tau <- seq_len(2 * B) / (2 * B)
  # Do any features exceed the upper bound?
  dat <- data.table(tau, err_bound = ub)[tau >= epsilon]
  dat[, detected := sapply(seq_len(nrow(dat)), function(i) sum(r >= dat$tau[i]))]
  dat[, surplus := ifelse(detected > err_bound, 1, 0)]
  # Export
  out <- data.table(
    'order' = order, 'rule' = rule, 
    'decision' = ifelse(sum(dat$surplus) > 0, 1, 0)
  )
  return(out)
}


