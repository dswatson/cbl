### Simulations for subgraph discovery with many ancestors ###

# Load libraries and Shah's code, register cores
source('shah_ss.R')
library(data.table)
library(pcalg)
library(RBGL)
library(matrixStats)
library(glmnet)
library(lightgbm)
library(tidyverse)
library(doMC)
registerDoMC(16)

# Set seed
set.seed(123, kind = "L'Ecuyer-CMRG")

################################################################################

### SIMULATION ###

#' @param n Sample size.
#' @param d_z Dimensionality of Z.
#' @param d_x Dimensionality of X.
#' @param rho Auto-correlation of the Toeplitz matrix for Z.
#' @param r2 Proportion of variance explained for all foreground variables. 
#' @param lin_pr Probability that an edge denotes a linear relationship.
#' @param method Method used for generating the graph structure. Options are
#'   \code{"er"} for Erdós-Rényi and \code{"barabasi"} for Barabási-Albert.
#' @param sp Average sparsity of the graph. Note that this must be high for 
#'   \code{method = "barabasi"} or else you'll run into errors.
#' @param pref Strength of preferential attachment if \code{method = "barabasi"}.
#' 

# Data simulation function
# Note: lower triangular adj_mat means that column is a parent of row
sim_dat <- function(n, d_z, d_x, rho, r2, lin_pr, sp, method, pref) {
  # Simulate background variables
  z <- matrix(rnorm(n * d_z), ncol = d_z)
  var_z <- 1 / d_z
  Sigma <- toeplitz(rho^(0:(d_z - 1))) * var_z
  z <- z %*% chol(Sigma)
  colnames(z) <- paste0('z', seq_len(d_z))
  # Optionally apply nonlinear transformations
  prep <- function(dat, pr) {
     out <- dat
     if (pr < 1) {
       # Pick features to transform
       n_nl <- round((1 - pr) * ncol(dat))
       if (n_nl > 0) {
         tmp <- data.table(idx = sample.int(ncol(dat), size = n_nl))
         tmp[, nl := sample(c('sq', 'sqrt', 'sftpls', 'relu'), 
                            size = n_nl, replace = TRUE)]
         out[, tmp[nl == 'sq', idx]] <- dat[, tmp[nl == 'sq', idx]]^2
         out[, tmp[nl == 'sqrt', idx]] <- sqrt(abs(dat[, tmp[nl == 'sqrt', idx]]))
         out[, tmp[nl == 'sftpls', idx]] <- log(1 + exp(dat[, tmp[nl == 'sftpls', idx]]))
         out[, tmp[nl == 'relu', idx]] <- ifelse(dat[, tmp[nl == 'relu', idx]] > 0, 
                                                 dat[, tmp[nl == 'relu', idx]], 0)
       }
     }
     return(out)
  }
  # Generate noise
  sim_noise <- function(signal, r2) {
    var_mu <- var(signal)
    var_noise <- (var_mu - r2 * var_mu) / r2
    noise <- rnorm(n, sd = sqrt(var_noise))
    return(noise)
  }
  # Simulate graph
  m <- (1 - sp) * (d_z + d_x - 1)
  g <- randDAG(d_z + d_x, m, method = method, par1 = pref, weighted = FALSE)
  t_srt <- as.numeric(tsort(g))
  z_idx <- data.table(z = seq_len(d_z), g = t_srt[seq_len(d_z)])
  x_idx <- data.table(x = seq_len(d_x), g = t_srt[(d_z + 1):(d_z + d_x)])
  # Compute X recursively, record adjacency matrix
  x_labs <- paste0('x', seq_len(d_x))
  x <- matrix(nrow = n, ncol = d_x, dimnames = list(NULL, x_labs))
  adj_mat <- matrix(0, nrow = d_x, ncol = d_x, dimnames = list(x_labs, x_labs))
  diag(adj_mat) <- NA_real_
  for (j in seq_len(d_x)) {
    # Index the parents
    pa <- c()
    for (i in seq_len(d_z + j - 1)) {
      if (any(grepl(t_srt[d_z + j], as.numeric(g@edgeL[[t_srt[i]]]$edges)))) {
        pa <- c(pa, t_srt[i])
      }
    }
    # Compute Z signal with Rademacher weights
    pa_z <- prep(z[, z_idx[g %in% pa, z]], lin_pr)
    beta_z <- sample(c(1, -1), size = ncol(pa_z), replace = TRUE)
    signal_z <- as.numeric(pa_z %*% beta_z)
    # Compute X signal, if applicable
    if (any(x_idx$g %in% pa)) {
      pa_x <- as.matrix(prep(x[, x_idx[g %in% pa, x]], lin_pr))
      adj_mat[j, x_idx[g %in% pa, x]] <- 1
      causal_wt <- 1 / length(pa)
      sigma_xij <- sqrt(causal_wt * var(signal_z))
      beta_x <- sigma_xij / colSds(pa_x)
      signal_x <- as.numeric(pa_x %*% beta_x)
    } else {
      signal_x <- 0
    }
    signal_xj <- signal_z + signal_x
    # Add appropriate noise and export
    x[, j] <- signal_xj + sim_noise(signal_xj, r2)
  }
  # Export
  params <- list(
    'n' = n, 'd_z' = d_z, 'd_x' = d_x, 'rho' = rho, 'r2' = r2, 'lin_pr' = lin_pr, 
    'sp' = sp, 'method' = method, 'pref' = pref
  )
  out <- list('dat' = data.table(z, x), 'adj_mat' = adj_mat, 'params' = params)
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


#' @param df Table of (de)activation rates.
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


#' @param df Table of (de)activation rates.
#' @param eps Consistency lower bound, as computed by \code{epsilon_fn}.
#' @param order Causal order of interest, either \code{"ij"} or \code{"ji"}.
#' @param rule Inference rule, either \code{"R1"} or \code{"R2"}.
#' @param B Number of complementary pairs to draw for stability selection.

# Infer causal direction using stability selection
ss_fn <- function(df, eps, order, rule, B) {
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
  dat <- data.frame(tau, err_bound = ub) %>%
    filter(tau >= eps) %>%
    rowwise() %>%
    mutate(detected = sum(r >= tau)) %>% 
    ungroup(.) %>%
    mutate(surplus = ifelse(detected > err_bound, 1, 0))
  # Export
  out <- data.table(
    'order' = order, 'rule' = rule, 
    'decision' = ifelse(sum(dat$surplus) > 0, 1, 0)
  )
  return(out)
}

#' @param b Subsample index.
#' @param i First foreground variable index.
#' @param j Second foreground variable index.
#' @param x Matrix of foreground variables.
#' @param z_t Intersection of iteration-t known non-descendants for foreground
#'   variables \code{i} and \code{j}.
#' @param s Regression method. 
#' @param params Optional list of parameters to use when \code{s = "boost"}.

# Complementary pairs subsampling loop
sub_loop <- function(b, i, j, x, z_t, s, params) {
  # Prelimz
  n <- nrow(x) 
  d_zt <- ncol(z_t)
  # Take complementary subsets
  a_set <- sample(n, round(0.5 * n))
  b_set <- seq_len(n)[-a_set]
  # Fit reduced models
  s0 <- sapply(c(i, j), function(k) {
    c(l0(z_t[a_set, ], x[a_set, k], s, params), 
      l0(z_t[b_set, ], x[b_set, k], s, params))
  })
  # Fit expanded models
  s1 <- sapply(c(i, j), function(k) {
    not_k <- setdiff(c(i, j), k)
    c(l0(cbind(z_t[a_set, ], x[a_set, not_k]), x[a_set, k], s, params),
      l0(cbind(z_t[b_set, ], x[b_set, not_k]), x[b_set, k], s, params))
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


#' @param sim_obj Simulation object as computed by \code{sim_dat}.
#' @param maxiter Maximum number of iterations to loop through if convergence
#'   is elusive.
#' @param gamma Omission threshold.
#' @param B Number of complementary pairs to draw for stability selection.

# Subdag discovery via confounder blanket regression
cbl_fn <- function(sim_obj, gamma = 0.5, maxiter = 100, B = 50) {
  ### PRELIMINARIES ###
  # Get data, hyperparameters, train/test split
  dat <- sim_obj$dat
  n <- nrow(dat)
  z <- as.matrix(select(dat, starts_with('z')))
  d_z <- ncol(z)
  x <- as.matrix(select(dat, starts_with('x')))
  d_x <- ncol(x)
  xlabs <- paste0('x', seq_len(d_x))
  linear <- sim_obj$params$lin_pr == 1
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
            df <- foreach(bb = seq_len(B), .combine = rbind) %do%
              sub_loop(bb, i, j, x, z_t, f, prms)
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
              if (out[rule == 'R1' & order == 'ji', decision == 1]) {
                adj1[i, j] <- 1
                adj1[j, i] <- 0
              } else if (out[rule == 'R1' & order == 'ij', decision == 1]) {
                adj1[j, i] <- 1
                adj1[i, j] <- 0
              } else if (out[rule == 'R2', sum(decision) == 2]) {
                adj1[i, j] <- adj1[j, i] <- 0
              } else if (out[rule == 'R2' & order == 'ji', decision == 1]) {
                adj1[i, j] <- 0.5 
              } else if (out[rule == 'R2' & order == 'ij', decision == 1]) {
                adj1[j, i] <- 0.5
              } 
            }
          }
        } 
      }
    }
    # Check for transitivity
    closure <- FALSE
    while(closure == FALSE) {
      closure <- TRUE
      for (i in seq_len(d_x)) {
        m <- which(adj1[, i] == 1)
        e <- which(adj1[, m] == 1)
        if (length(e) > 0) {
          if (is.na(adj1[e, i]) | adj1[e, i] != 1) {
            adj1[e, i] <- 1
            closure <- FALSE
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

################################################################################

### RFCI ###

rfci_fn <- function(sim_obj) {
  # Extract data
  dat <- sim_obj$dat
  n <- nrow(dat)
  d_z <- sim_obj$params$d_z
  zlabs <- paste0('z', seq_len(d_z))
  d_x <- sim_obj$params$d_x
  xlabs <- paste0('x', seq_len(d_x))
  k <- round(sim_obj$params$sp * d_z)
  # Gap matrix ensures we don't compute intra-Z edges
  rng <- (d_z + 1):(d_z + d_x)
  gps <- matrix(TRUE, nrow = d_z + d_x, ncol = d_z + d_x)
  gps[rng, ] <- gps[, rng] <- FALSE
  # RFCI (returns PAG)
  rho_list <- list(C = cor(dat), n = n)
  # This is not the original (skel.method = 'original'), but it's the only one 
  # that completes in a reasonable amount of time
  rfci_out <- rfci(rho_list, indepTest = gaussCItest, alpha = 0.1, 
                   labels = c(zlabs, xlabs), 
                   skel.method = 'stable.fast', numCores = 8,
                   fixedGaps = gps, m.max = k)
  rfci_amat <- rfci_out@amat[rng, rng]
  # Export
  return(rfci_amat)
}

################################################################################

### GES ###

ges_fn <- function(sim_obj) {
  # Extract data
  dat <- sim_obj$dat
  n <- nrow(dat)
  d_z <- sim_obj$params$d_z
  d_x <- sim_obj$params$d_x
  xlabs <- paste0('x', seq_len(d_x))
  k <- round(sim_obj$params$sp * d_z)
  # Gap matrix ensures we don't compute intra-Z edges
  rng <- (d_z + 1):(d_z + d_x)
  gps <- matrix(TRUE, nrow = d_z + d_x, ncol = d_z + d_x)
  gps[rng, ] <- gps[, rng] <- FALSE
  # GES (returns CPDAG)
  score <- new('GaussL0penObsScore', dat) # BCI
  ges_out <- ges(score, labels = score$getNodes(), maxDegree = k,
                 fixedGaps = gps, phase = c('forward', 'backward'),
                 iterate = FALSE) # Original Chickering algorithm
  in_edges <- ges_out$essgraph$.in.edges[rng]
  in_edges <- lapply(seq_along(in_edges), function(k) in_edges[[k]] - d_z)
  ges_amat <- matrix(0, nrow = d_x, ncol = d_x, dimnames = list(xlabs, xlabs))
  for (i in seq_len(d_x)) {
    for (j in seq_len(d_x)[-i]) {
      if (j %in% in_edges[[i]]) {
        ges_amat[i, j] <- 1
      }
    }
  }
  return(ges_amat)
}


################################################################################

# Initialize simulations
out <- data.table(
  n = NA, d_z = NA, sp = NA, idx = NA, amat_cbl = NA, amat_ges = NA, amat = NA
)
saveRDS(out, './res/multivar.rds')

big_loop <- function(n, d_z, sp, i) {
  # Simulate data
  sim_obj <- sim_dat(n = n, d_z = d_z, d_x = 6, rho = 1/4, 
                     r2 = 2/3, lin_pr = 1, 
                     sp = sp, method = 'er', pref = 1)
  # Estimate ancestor matrix via CBL
  amat_cbl <- cbl_fn(sim_obj)
  # Estimate CPDAG via GES
  amat_ges <- ges_fn(sim_obj)
  # Export results
  new <- data.table(
    'n' = n, 'd_z' = d_z, 'sp' = sp, 'idx' = i, 
    amat_cbl = list(amat_cbl), 
    amat_ges = list(amat_ges),
    amat = list(sim_obj$adj_mat)
  )
  old <- readRDS('./res/multivar.rds')
  out <- na.omit(rbind(old, new))
  saveRDS(out, './res/multivar.rds')
}

# Compute in parallel
foreach(nn = c(500, 1000, 2000, 4000, 8000)) %:%
  foreach(dd = c(50, 100)) %:%
  foreach(ss = c(1/4, 3/4)) %:%
  foreach(ii = 1:20) %dopar%
  big_loop(nn, dd, ss, ii)


# RFCI loop computed separately
out <- data.table(
  n = NA, d_z = NA, sp = NA, idx = NA, amat_rfci = NA, amat = NA
)
saveRDS(out, './res/multivar_rfci.rds')
rfci_loop <- function(n, d_z, sp, i) {
  # Simulate data
  sim_obj <- sim_dat(n = n, d_z = d_z, d_x = 6, rho = 1/4, 
                     r2 = 2/3, lin_pr = 1, 
                     sp = sp, method = 'er', pref = 1)
  # Run RFCI
  amat_rfci <- rfci_fn(sim_obj)
  # Export results
  new <- data.table(
    'n' = n, 'd_z' = d_z, 'sp' = sp, 'idx' = i,
    amat_rfci = list(amat_rfci), 
    amat = list(sim_obj$adj_mat)
  )
  old <- readRDS('./res/multivar_rfci.rds')
  out <- na.omit(rbind(old, new))
  saveRDS(out, './res/multivar_rfci.rds')
}
foreach(dd = c(50, 100)) %:%
  foreach(ss = c(1/4, 3/4)) %:%
  foreach(ii = 1:5) %do%
  rfci_loop(n = 500, dd, ss, ii)
foreach(ss = c(1/4, 3/4)) %:%
  foreach(ii = 1:5) %do%
  rfci_loop(n = 1000, d_z = 50, ss, ii)





