# Load libraries and Shah's code, register cores
source('shah_ss.R')
library(data.table)
library(trigger)
library(glmnet)
library(bestsubset)
library(tidyverse)
library(doMC)
registerDoMC(8)

# Set seed
set.seed(123, kind = "L'Ecuyer-CMRG")

#' @param x Design matrix.
#' @param y Outcome vector.
#' @param trn Training indices.
#' @param tst Test indices.
#' @param f Regression method.

# Fit regressions, return bit vector for feature selection.
l0 <- function(x, y, trn, tst, f) {
  out <- double(length = ncol(x))
  names(out) <- colnames(x)
  x_trn <- as.data.frame(x[trn, ])
  x_trn <- as.matrix(x_trn[!duplicated(as.list(x_trn))])
  x_tst <- x[tst, colnames(x_trn)]
  if (f == 'lasso') {
    fit <- glmnet(x_trn, y[trn], intercept = FALSE)
    y_hat <- predict(fit, newx = x_tst, s = fit$lambda)
    betas <- coef(fit, s = fit$lambda)[-1, ]
  } else if (f == 'step') {
    fit <- fs(x_trn, y[trn], intercept = FALSE, verbose = FALSE)
    y_hat <- predict(fit, newx = x_tst)
    betas <- coef(fit)[1:ncol(x_trn), ]
  } 
  epsilon <- y_hat - y[tst]
  mse <- colMeans(epsilon^2)
  betas <- betas[, which.min(mse)]
  nonzero <- names(betas)[betas != 0]
  out[nonzero] <- 1
  return(out)
}


#' @param df Table of (de)activation rates.
#' @param B Number of complementary pairs to draw for stability selection.

# Compute consistent lower bound
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
#' @param lb Consistency lower bound, as computed by \code{lb_fn}.
#' @param order Causal order of interest, either \code{"ij"} or \code{"ji"}.
#' @param rule Inference rule, either \code{"R1"} or \code{"R2"}.
#' @param B Number of complementary pairs to draw for stability selection.

# Infer causal direction using stability selection
ss_fn <- function(df, lb, order, rule, B) {
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
    filter(tau >= lb) %>%
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


#' @param g Target gene.
#' @param bp Size of window in base pairs.

# Return candidate eQTLs for a given gene and window
eqtl_cands <- function(g, bp = 5000) {
  chrom <- x.pos[gene == g, chr]
  lb <- x.pos[gene == g, lo] - bp
  ub <- x.pos[gene == g, hi] + bp
  out <- z.pos[chr == chrom & pos >= lb & pos <= ub, idx]
  return(out)
}


#' @param z Matrix of genotype data.
#' @param x Matrix of gene expression data.
#' @param z.pos Position of SNPs.
#' @param x.pos Position of genes.
#' @param f Regression method.
#' @param gamma Omission threshold.
#' @param maxiter Maximum number of iterations.
#' @param B Number of complementary pairs to draw for stability selection.

# Subdag discovery algorithm
subdag <- function(z, x, z.pos, x.pos, f, gamma = 0.5, maxiter = Inf, B = 50) {
  ### PRELIMINARIES ###
  # Get data, hyperparameters, train/test split
  n <- nrow(z)
  d_x <- ncol(x)
  # Initialize
  adj_list <- list(
    matrix(NA_real_, nrow = d_x, ncol = d_x, 
           dimnames = list(colnames(x), colnames(x)))
  )
  converged <- FALSE
  iter <- 0
  ### LOOP IT ###
  while(converged == FALSE & iter <= maxiter) {
    # Extract relevant adjacency matrices
    if (iter == 0) {
      adj0 <- adj1 <- adj_list[[1]]
    } else {
      adj0 <- adj_list[[iter]]
      adj1 <- adj_list[[iter + 1]]
    }
    # Subsampling loop
    sub_loop <- function(b, i, j, a_t) {
      d_at <- ncol(a_t)
      # Take complementary subsets
      a_set <- sample(n, round(0.5 * n))
      a_trn <- sample(a_set, round(0.8 * length(a_set)))
      a_tst <- setdiff(a_set, a_trn)
      b_set <- seq_len(n)[-a_set]
      b_trn <- sample(b_set, round(0.8 * length(b_set)))
      b_tst <- setdiff(b_set, b_trn)
      # Fit reduced models
      s0 <- sapply(c(i, j), function(k) {
        c(l0(a_t, x[, k], a_trn, a_tst, f), 
          l0(a_t, x[, k], b_trn, b_tst, f))
      })
      # Fit expanded models
      s1 <- sapply(c(i, j), function(k) {
        not_k <- setdiff(c(i, j), k)
        a_tmp <- cbind(a_t, x[, not_k])
        colnames(a_tmp)[ncol(a_tmp)] <- colnames(x)[not_k]
        c(l0(a_tmp, x[, k], a_trn, a_tst, f),
          l0(a_tmp, x[, k], b_trn, b_tst, f))
      })
      # Record disconnections and (de)activations
      dis_a <- any(s1[d_at + 1, ] == 0)
      dis_b <- any(s1[2 * (d_at + 1), ] == 0)
      dis <- rep(c(dis_a, dis_b), each = d_at)
      d_ji <- s0[, 1] == 1 & s1[seq_len(d_at), 1] == 0
      a_ji <- s0[, 2] == 0 & s1[seq_len(d_at), 2] == 1
      d_ij <- s0[, 2] == 1 & s1[seq_len(d_at), 2] == 0
      a_ij <- s0[, 1] == 0 & s1[seq_len(d_at), 1] == 1
      # Export
      out <- data.table(b = rep(c(2 * b - 1, 2 * b), each = d_at), i, j,
                        z = rep(colnames(a_t), times = 2),
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
          # Only continue if the set of nondescendants has increased since last 
          # iteration (i.e., have we learned anything new?)
          if (iter == 0 | length(a1) > length(a0)) {
            z_i <- eqtl_cands(g = colnames(x)[i])
            z_j <- eqtl_cands(g = colnames(x)[j])
            #z_t <- cbind(z[, intersect(z_i, z_j)], x[, a1])
            a_t <- cbind(z[, c(z_i, z_j)], x[, a1])
            df <- foreach(bb = seq_len(B), .combine = rbind) %dopar%
              sub_loop(bb, i, j, a_t)
            #print(c(i, j, iter))
            
            # Compute rates
            df[, disr := sum(dis) / .N]
            if (df$disr[1] >= gamma) { 
              adj1[i, j] <- adj1[j, i] <- 0
            } else {
              df[, drji := sum(d_ji) / .N, by = z]
              df[, arji := sum(a_ji) / .N, by = z]
              df[, drij := sum(d_ij) / .N, by = z]
              df[, arij := sum(a_ij) / .N, by = z]
              df <- unique(df[, .(i, j, z, disr, drji, arji, drij, arij)])
              # Consistent lower bound
              lb <- epsilon_fn(df, B)
              # Stable upper bound
              out <- foreach(oo = c('ji', 'ij'), .combine = rbind) %:%
                foreach(rr = c('R1', 'R2'), .combine = rbind) %do%
                ss_fn(df, lb, oo, rr, B)
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

################################################################################

# Import data
data(yeast)

# Extract, format data
z <- t(yeast$marker - 1)
colnames(z) <- paste0('z', 1:ncol(z))
x <- t(yeast$exp)
genes <- colnames(x)
x.pos <- as.data.table(yeast$exp.pos)
colnames(x.pos) <- c('chr', 'lo', 'hi')
x.pos[, gene := genes]
z.pos <- as.data.table(yeast$marker.pos)
colnames(z.pos) <- c('chr', 'pos')
z.pos[, idx := .I]

# Weirdly, there's some inconsistent naming in the data, so YJR008W = MHO1
sub_g <- c('CHO1', 'ITR1', 'YJR008W', 'OPI3', 'OPT1', 'THI7')
x_try <- x[, sub_g]
a_mat <- subdag(z, x_try, z.pos, x.pos, 
                f = 'lasso', gamma = 0.5, maxiter = Inf, B = 50)

# Plot results
library(igraph)
library(ggsci)
edge_df <- data.frame(
  from = c('ITR1', 'ITR1', 'MHO1', 'MHO1', 'THI7', 'THI7', 'OPT1', 'OPI3'),
  to = c('OPI3', 'CHO1', 'OPI3', 'CHO1', 'OPT1', 'CHO1', 'CHO1', 'CHO1')
)
g <- graph_from_data_frame(edge_df)
V(g)$label.cex = 2
V(g)$label.color = 'black'
plot(g, vertex.size = 50, vertex.color = pal_jco()(5)[5], 
     edge.color = 'black', layout = layout_nicely)


