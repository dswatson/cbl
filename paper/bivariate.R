### Simulations for subgraph discovery with many ancestors ###

# Load libraries and Shah's code, register cores
source('shah_ss.R')
library(data.table)
library(glmnet)
library(lightgbm)
library(ppcor)
library(tidyverse)
library(doMC)
registerDoMC(16)

# Set seed
set.seed(123, kind = "L'Ecuyer-CMRG")

################################################################################

### SIMULATION ###

#' @param n Sample size.
#' @param d_z Dimensionality of Z.
#' @param rho Auto-correlation of the Toeplitz matrix for Z.
#' @param sp Sparsity of the connections from background to foreground.
#' @param r2 Proportion of variance explained by endogenous features.
#' @param lin_pr Probability that an edge denotes a linear relationship.
#' @param g Expected output of an independence oracle, one of \code{"xy"},
#'   \code{"ci"}, or \code{"na"}.
#' 

# Data simulation function
sim_dat <- function(n, d_z, rho, sp, r2, lin_pr, g) {
  # Simulate ancestors
  z <- matrix(rnorm(n * d_z), ncol = d_z)
  var_z <- 1 / d_z # Does this make a difference?
  Sigma <- toeplitz(rho^(0:(d_z - 1))) * var_z
  z <- z %*% chol(Sigma)
  colnames(z) <- paste0('z', seq_len(d_z))
  # Random Rademacher weights
  beta <- gamma <- double(length = d_z)
  k <- round((1 - sp) * d_z)
  beta[sample(d_z, k)] <- sample(c(1, -1), k, replace = TRUE)
  gamma[sample(d_z, k)] <- sample(c(1, -1), k, replace = TRUE)
  # Nonlinear transformations?
  if (lin_pr < 1) {
    z_out <- zz <- z
    # Create matrix zz of random nonlinear transformations
    idx <- split(sample.int(d_z), sort(seq_len(d_z) %% 4))
    names(idx) <- c('sq', 'sqrt', 'sftpls', 'relu')
    zz[, idx$sq] <- z[, idx$sq]^2
    zz[, idx$sqrt] <- sqrt(abs(z[, idx$sqrt]))
    zz[, idx$sftpls] <- log(1 + exp(z[, idx$sftpls]))
    zz[, idx$relu] <- ifelse(z[, idx$relu] > 0, z[, idx$relu], 0)
    # Sample columns from zz with probability 1 - lin_pr
    nonlin_idx <- sample.int(d_z, d_z * (1 - lin_pr))
    z_out[, nonlin_idx] <- zz[, nonlin_idx]
    signal_x <- as.numeric(z_out %*% beta)
    signal_y <- as.numeric(z_out %*% gamma)
  } else {
    signal_x <- as.numeric(z %*% beta)
    signal_y <- as.numeric(z %*% gamma)
  }
  # Generate noise
  sim_noise <- function(signal, r2) {
    var_mu <- var(signal)
    var_noise <- (var_mu - r2 * var_mu) / r2
    noise <- rnorm(n, sd = sqrt(var_noise))
    return(noise)
  }
  # X data
  x <- signal_x + sim_noise(signal_x, r2)
  # Identifiable?
  if (g == 'na') {
    shared_parents <- which(beta != 0 & gamma != 0)
    u_idx <- sample(shared_parents, size = length(shared_parents)/2)
    z <- z[, -u_idx]
    d_z <- ncol(z)
    d_u <- length(u_idx)
    beta <- beta[-u_idx]
    gamma <- gamma[-u_idx]
  } else {
    d_u <- 0
  }
  # Y data
  if (g %in% c('ci', 'na')) {
    y <- signal_y + sim_noise(signal_y, r2)
  } else if (g == 'xy') {
    signal_z_to_y <- signal_y
    xzr <- 1 / (k + 1)
    sigma_xy <- sqrt(xzr * var(signal_z_to_y))
    gamma_x <- sigma_xy / sd(x)
    gamma <- c(gamma, gamma_x)
    signal_y <- signal_z_to_y + x * gamma_x
    y <- signal_y + sim_noise(signal_y, r2)
  }
  # Export
  params <- list(
    'n' = n, 'd_z' = d_z, 'd_u' = d_u, 'rho' = rho, 'sp' = sp, 'r2' = r2, 
    'lin_pr' = lin_pr, 'g' = g
  )
  out <- list(
    'dat' = data.table(z, 'x' = x, 'y' = y),
    'wts' = list('beta' = beta, 'gamma' = gamma), 'params' = params
  )
  return(out)
}

################################################################################

### CONFOUNDER BLANKET LEARNER ###

#' @param x Design matrix.
#' @param y Outcome vector.
#' @param f Regression method, either \code{"lasso"} or \code{"gbm"}.
#' @param prms List of parameters to use when \code{f = "gbm"}.
#' 

# Fit regressions, return bit vector for feature selection.
l0 <- function(x, y, f, prms) {
  n <- nrow(x)
  trn <- sample(n, round(0.8 * n))
  tst <- seq_len(n)[-trn]
  if (f == 'lasso') {
    fit <- glmnet(x[trn, ], y[trn], intercept = FALSE)
    y_hat <- predict(fit, newx = x[tst, ], s = fit$lambda)
    eps <- y_hat - y[tst]
    mse <- colMeans(eps^2)
    betas <- coef(fit, s = fit$lambda)[-1, which.min(mse)]
    out <- ifelse(betas == 0, 0, 1)
  } else if (f == 'gbm') {
    d_trn <- lgb.Dataset(x[trn, ], label = y[trn])
    d_tst <- lgb.Dataset.create.valid(d_trn, x[tst, ], label = y[tst])
    fit <- lgb.train(params = prms, data = d_trn, valids = list(tst = d_tst), 
                     nrounds = 3500, early_stopping_rounds = 10, verbose = 0)
    vimp <- lgb.importance(fit)
    out <- as.numeric(colnames(x) %in% vimp$Feature)
  } 
  return(out)
}


#' @param z Matrix of background variables.
#' @param x Candidate cause.
#' @param y Candidate effect.
#' @param linear Are all structural equations linear?
#' @param B Number of complementary pairs to draw for stability selection.
#' 

# Compute (de)activation rates for X -> Y and Y -> X
rate_fn <- function(z, x, y, linear, B) {
  # Preliminaries
  n <- nrow(z)
  d_z <- ncol(z)
  zx <- cbind(z, x)
  zy <- cbind(z, y)
  if (linear) {
    f <- 'lasso'
    prms <- NULL
  } else {
    f <- 'gbm'
    prms <- list(
      objective = 'regression', max_depth = 1, 
      bagging.fraction = 0.5, feature_fraction = 0.8, 
      num_threads = 1, force_col_wise = TRUE
    )
  }
  # Compute disconnections and (de)activations per subsample
  fit_fn <- function(b) {
    # Take complementary subsets
    i_set <- sample(n, round(0.5 * n))
    j_set <- seq_len(n)[-i_set]
    # Compute active sets
    s <- data.frame(
      y0 = c(l0(z[i_set, ], y[i_set], f, prms), NA_real_, 
             l0(z[j_set, ], y[j_set], f, prms), NA_real_), 
      y1 = c(l0(zx[i_set, ], y[i_set], f, prms), 
             l0(zx[j_set, ], y[j_set], f, prms)),
      x0 = c(l0(z[i_set, ], x[i_set], f, prms), NA_real_, 
             l0(z[j_set, ], x[j_set], f, prms), NA_real_), 
      x1 = c(l0(zy[i_set, ], x[i_set], f, prms), 
             l0(zy[j_set, ], x[j_set], f, prms))
    )
    # Record disconnections and (de)activations
    dis_i <- any(c(s$y1[d_z + 1], s$x1[d_z + 1]) == 0)
    dis_j <- any(c(s$y1[2 * (d_z + 1)], s$x1[2 * (d_z + 1)]) == 0)
    dis <- rep(c(dis_i, dis_j), each = d_z + 1)
    d_xy <- s$y0 == 1 & s$y1 == 0
    a_xy <- s$x0 == 0 & s$x1 == 1
    d_yx <- s$x0 == 1 & s$x1 == 0
    a_yx <- s$y0 == 0 & s$y1 == 1
    extras <- c(d_z + 1, 2 * (d_z + 1))
    d_xy[extras] <- a_xy[extras] <- d_yx[extras] <- a_yx[extras] <- NA_real_
    # Export
    out <- data.table(b = rep(c(2 * b - 1, 2 * b), each = d_z + 1), 
                      dis, d_xy, a_xy, d_yx, a_yx,
                      z = rep(seq_len(d_z + 1), times = 2))
    return(out)
  }
  out <- foreach(bb = seq_len(B), .combine = rbind) %do% 
    fit_fn(bb)
  # Compute rates
  out[, disr := sum(dis) / .N]
  out[, drxy := sum(d_xy) / .N, by = z]
  out[, arxy := sum(a_xy) / .N, by = z]
  out[, dryx := sum(d_yx) / .N, by = z]
  out[, aryx := sum(a_yx) / .N, by = z]
  # Tidy up, export
  out <- unique(out[, .(z, disr, drxy, arxy, dryx, aryx)])
  return(out)
}


#' @param df Results object output by \code{rate_fn}.
#' @param B Number of complementary pairs to draw for stability selection.

# Compute consistency lower bound
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


#' @param res Results object output by \code{rate_fn}.
#' @param eps Lower bound output by \code{epsilon_fn}.
#' @param order Assume X \preceq Y or Y \preceq X?
#' @param rule Detect via deactivation (\code{"R1"}) or activation (\code{"R2"})?
#' @param B Number of complementary pairs to draw for stability selection.
#' 

# Infer causal direction using stability selection
ss_fn <- function(res, eps, order, rule, B) {
  # Subset the data
  if (order == 'xy' & rule == 'R1') {
    r <- res$drxy 
  } else if (order == 'xy' & rule == 'R2') {
    r <- res$arxy 
  } else if (order == 'yx' & rule == 'R1') {
    r <- res$dryx 
  } else if (order == 'yx' & rule == 'R2') {
    r <- res$aryx 
  }
  r <- na.omit(r)
  if (max(r) == 0) {
    dat <- data.frame(surplus = 0)
  } else {
    # Stability selection parameters
    theta <- mean(r)
    ub <- minD(theta, B) * sum(r <= theta)
    tau <- seq_len(2 * B) / (2 * B)
    # Do any features exceed the upper bound?
    dat <- data.frame(tau, err_bound = ub) %>%
      filter(tau >= eps) %>%
      rowwise() %>%
      mutate(detected = sum(r >= tau)) %>% 
      ungroup() %>%
      mutate(surplus = ifelse(detected > err_bound, 1, 0))
  }
  # Export
  out <- data.table(
    'order' = order, 'rule' = rule, 
    'decision' = ifelse(sum(dat$surplus) > 0, 1, 0)
  )
  return(out)
}


#' @param z Matrix of background variables.
#' @param x Candidate cause.
#' @param y Candidate effect.
#' @param linear Are all structural equations linear?
#' @param gamma Omission threshold.
#' 

# Wrap it up
cbl_fn <- function(z, x, y, linear, gamma) {
  # Compute rates for each z
  res <- rate_fn(z, x, y, linear, B = 50)
  # Disconnected?
  if (res$disr[1] > gamma) {
    decision <- 'ci'
  } else {
    # (De)activation rates
    eps <- epsilon_fn(res, B = 50)
    sum_tbl <- foreach(oo = c('xy', 'yx'), .combine = rbind) %:%
      foreach(rr = c('R1', 'R2'), .combine = rbind) %do%
      ss_fn(res, eps, oo, rr, B = 50)
    if (sum_tbl[rule == 'R2', sum(decision) == 2]) {
      decision <- 'ci'
    } else if (sum_tbl[order == 'xy', sum(decision) > 0]) {
      decision <- 'xy'
    } else if (sum_tbl[order == 'yx', sum(decision) > 0]) {
      decision <- 'yx'
    } else {
      decision <- 'na'
    }
  }
  # Export
  out <- data.table(method = 'cbl', g_hat = decision) 
  return(out)
}

################################################################################

### ENTNER METHOD ###

#' @param x First vector.
#' @param y Second vector.
#' @param z Conditioning set.
#' @param trn Training index.
#' @param val Validation index.
#' @param tst Test index.
#' @param prms List of model parameters.
#' 

# Generalized covariance measure (Shah & Peters, 2020). 
# Tests conditional independence of x and y given z using gradient boosting.
gcm_test <- function(x, y, z, trn, val, tst, prms) {
  # Model 1
  d1_trn <- lgb.Dataset(z[trn, ], label = y[trn])
  d1_val <- lgb.Dataset.create.valid(d1_trn, z[val, ], label = y[val])
  f1 <- lgb.train(params = prms, data = d1_trn, valids = list(val = d1_val),
                  nrounds = 3000, early_stopping_rounds = 10, verbose = 0)
  eps1 <- y[tst] - predict(f1, z[tst, ])
  # Model 2
  d2_trn <- lgb.Dataset(z[trn, ], label = x[trn])
  d2_val <- lgb.Dataset.create.valid(d2_trn, z[val, ], label = x[val])
  f2 <- lgb.train(params = prms, data = d2_trn, valids = list(val = d2_val),
                  nrounds = 3000, early_stopping_rounds = 10, verbose = 0)
  eps2 <- x[tst] - predict(f2, z[tst, ])
  # Inference
  nn <- length(tst)
  R <- eps1 * eps2
  R.sq <- R^2
  meanR <- mean(R)
  z_score <- sqrt(nn) * meanR / sqrt(mean(R.sq) - meanR^2)
  p_value <- 2 * pnorm(abs(z_score), lower.tail = FALSE)
  return(p_value)
}


#' @param z Matrix of background variables.
#' @param x Candidate cause.
#' @param y Candidate effect.
#' @param linear Are all structural equations linear?
#' @param k Size of expected admissible set.
#' @param alpha Significance threshold for inferring dependence.
#' @param tau Other threshold for "inferring" independence.
#' 

# Constraint-based: for this comparison, we presume X \preceq Y
# and assume access to the true data sparsity
constr_fn <- function(z, x, y, linear, k, alpha, tau) {
  # Preliminaries
  n <- nrow(z)
  d_z <- ncol(z)
  if (linear) {
    B <- 1000
    r1_thresh <- 1/200
  } else {
    B <- 500
    r1_thresh <- 0
    prms <- list(
      objective = 'regression', max_depth = 1, 
      bagging.fraction = 0.5, feature_fraction = 0.8, 
      num_threads = 1, force_col_wise = TRUE
    )
  }
  # Evaluate R1 10x more frequently than R2
  r2_idx <- seq_len(B/10)
  # Entner function
  entner <- function(b) {
    # Sample a random subset Z_b and variable W
    z_idx <- sample(d_z, k)
    z_b <- z[, z_idx]
    j <- sample(seq_len(d_z)[-z_idx], 1)
    w <- z[, j]
    if (linear) {
      ### RULE 1 ###
      pmat0 <- suppressWarnings(pcor(cbind(y, w, z_b))$p.value)
      p1.i <- pmat0[1, 2]
      pmat1 <- suppressWarnings(pcor(cbind(y, w, z_b, x))$p.value)
      p1.ii <- pmat1[1, 2]
      ### RULE 2 ###
      if (b %in% r2_idx) {
        pmat2 <- suppressWarnings(pcor(cbind(y, x, z_b))$p.value)
        p2.i <- pmat2[1, 2]
        pmat3 <- suppressWarnings(pcor(cbind(x, w, z_b))$p.value)
        p2.ii <- pmat3[1, 2]
      }
    } else {
      # Train/validation/test split
      trn <- sample(n, round(0.7 * n))
      val <- sample(setdiff(seq_len(n), trn), round(0.15 * n))
      tst <- seq_len(n)[-c(trn, val)]
      ### RULE 1 ###
      p1.i <- gcm_test(w, y, z_b, trn, val, tst, prms)
      p1.ii <- gcm_test(w, y, cbind(z_b, x), trn, val, tst, prms)
      ### RULE 2 ###
      if (b %in% r2_idx) {
        p2.i <- gcm_test(x, y, z_b, trn, val, tst, prms)
        p2.ii <- gcm_test(w, x, z_b, trn, val, tst, prms)
      } 
    }
    # Apply rules
    if (!b %in% r2_idx) {
      p2.i <- 0
      p2.ii <- 1
    }
    r1 <- ifelse(p1.i <= alpha & p1.ii >= tau, 1, 0)
    r2 <- ifelse(p2.i >= tau | (p2.ii <= alpha & p1.i >= tau), 1, 0)
    # Export
    out <- data.table(r1, r2)
    return(out)
  }
  # Apply Entner's rules with B random subset-variable pairs
  df <- foreach(bb = seq_len(B), .combine = rbind) %do%
    entner(bb)
  # Selecting different decision thresholds based on experimentation
  # Note -- this is very generous!
  if (df[, sum(r1)] > r1_thresh) {
    g_hat <- 'xy' 
  } else if (df[seq_len(B/10), sum(r2) / .N] >= 1/5) {
    g_hat <- 'ci'
  } else {
    g_hat <- 'na'
  }
  out <- data.table(method = 'constr', g_hat)
  return(out)
}

################################################################################

### SCORE-BASED METHOD ###

#' @param x Design matrix.
#' @param y Response vector.
#' @param trn Training index.
#' @param val Validation index.
#' @param tst Test index.
#' @param f Function class.
#' @param prms List of parameters to use when \code{f = "gbm"}.
#' 

# Perform feature selection on a training set, return a list with 
# (i) R^2 and (ii) test residuals
scr_fn <- function(x, y, trn, val, tst, f, prms) {
  # Fit models
  if (f == 'lasso') {
    fit <- glmnet(x[trn, ], y[trn])
    y_val <- predict(fit, newx = x[val, ], s = fit$lambda)
    eps <- y[val] - y_val
    mse <- colMeans(eps^2)
    y_hat <- predict(fit, newx = x[tst, ], s = fit$lambda[which.min(mse)])
  } else if (f == 'gbm') {
    d_trn <- lgb.Dataset(x[trn, ], label = y[trn])
    d_val <- lgb.Dataset.create.valid(d_trn, x[val, ], label = y[val])
    fit <- lgb.train(params = prms, data = d_trn, valids = list(val = d_val), 
                     nrounds = 3500, early_stopping_rounds = 10, verbose = 0)
    y_hat <- predict(fit, x[tst, ])
  }
  # Compute residuals, r2
  eps <- y[tst] - y_hat
  r2 <- 1 - (sum(eps^2) / sum((y[tst] - mean(y[tst]))^2))
  return(list('eps' = eps, 'r2' = r2))
}


#' @param z Matrix of background variables.
#' @param x Candidate cause.
#' @param y Candidate effect.
#' @param linear Are all structural equations linear?
#' @param alpha Significance threshold for testing.
#' 

# Evaluate graph structures using score-based method
score_fn <- function(z, x, y, linear, alpha) {
  # Preliminaries
  n <- nrow(z)
  d_z <- ncol(z)
  if (linear) {
    f <- 'lasso'
    prms <- NULL
  } else {
    f <- 'gbm'
    prms <- list(
      objective = 'regression', max_depth = 1, 
      bagging.fraction = 0.5, feature_fraction = 0.8, 
      num_threads = 1, force_col_wise = TRUE
    )
  }
  # Train/validation/test split
  trn <- sample(n, round(0.7 * n))
  val <- sample(setdiff(seq_len(n), trn), round(0.15 * n))
  tst <- seq_len(n)[-c(trn, val)]
  # Score X -> Y
  scr_xy <- scr_fn(cbind(z, x), y, trn, val, tst, f, prms)
  # Score Y -> X
  scr_yx <- scr_fn(cbind(z, y), x, trn, val, tst, f, prms)
  # Score X ~ Y
  scr_x <- scr_fn(z, x, trn, val, tst, f, prms)
  scr_y <- scr_fn(z, y, trn, val, tst, f, prms)
  # Decision procedure
  xy_scr <- scr_x$r2 + scr_xy$r2
  yx_scr <- scr_y$r2 + scr_yx$r2
  ci_scr <- scr_x$r2 + scr_y$r2
  if (ci_scr == max(c(xy_scr, yx_scr, ci_scr))) {
    g_hat <- 'ci'
  } else {
    if (xy_scr > yx_scr) {
      p_value <- cor.test(x[tst], scr_xy$eps)$p.value
      g_hat <- ifelse(p_value <= alpha, 'na', 'xy')
    } else {
      p_value <- cor.test(y[tst], scr_yx$eps)$p.value
      g_hat <- ifelse(p_value <= alpha, 'na', 'yx')
    }
  }
  # Export
  out <- data.table(method = 'score', g_hat)
  return(out)
}

################################################################################

### PUT IT ALL TOGETHER ###

# Initialize
out <- data.table(
  method = NA, g_hat = NA, n = NA, g = NA, idx = NA
)
saveRDS(out, './res/lin_biv_benchmark.rds')
saveRDS(out, './res/nl_biv_benchmark.rds')

# Big ol' wrapper
big_loop <- function(linear, n, g, i) {
  if (linear) {
    l_pr <- 1
    res_file <- './res/lin_biv_benchmark.rds'
  } else {
    l_pr <- 1/5
    res_file <- './res/nl_biv_benchmark.rds'
  }
  # Simulate data
  sim_obj <- sim_dat(n = n, d_z = 100, rho = 1/4, sp = 1/2, 
                     r2 = 2/3, lin_pr = l_pr, g = g)
  # Extract data
  dat <- sim_obj$dat
  z <- as.matrix(select(dat, starts_with('z')))
  d_z <- ncol(z)
  x <- dat$x
  y <- dat$y
  k <- round(d_z/4)
  # Confounder blanket learner
  df_b <- cbl_fn(z, x, y, linear, gamma = 0.5) 
  # Constraint function
  df_c <- constr_fn(z, x, y, linear, k, alpha = 0.1, tau = 0.5)
  # Score function
  df_s <- score_fn(z, x, y, linear, alpha = 0.1)
  # Import, export
  old <- readRDS(res_file)
  new <- rbind(df_b, df_c, df_s) %>%
    mutate('n' = n, 'g' = g, 'idx' = i) %>%
    as.data.table(.)
  out <- na.omit(rbind(old, new))
  saveRDS(out, res_file)
}

# Linear
foreach(ii = seq_len(100)) %:%
  foreach(nn = c(2500, 5000, 1e4)) %:%
  foreach(gg = c('xy', 'ci', 'na')) %dopar%
  big_loop(linear = TRUE, nn, gg, ii)

# Nonlinear
foreach(ii = seq_len(100)) %:%
  foreach(nn = c(5000, 1e4, 2e4)) %:%
  foreach(gg = c('xy', 'ci', 'na')) %dopar%
  big_loop(linear = FALSE, nn, gg, ii)







