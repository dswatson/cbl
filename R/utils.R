#' Feature selection subroutine
#' 
#' This function fits a potentially sparse supervised learning model and returns
#' a bit vector indicating which features were selected.
#' 
#' @param x Design matrix.
#' @param y Outcome vector.
#' @param s Regression method. Current options are \code{"lasso"} or
#'   \code{"boost"}.
#' @param params Optional list of parameters to use when \code{f = "boost"}.
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
      y_hat <- predict(fit, newx = x[tst, ], s = fit$lambda)
      eps <- y_hat - y[tst]
      mse <- colMeans(eps^2)
      betas <- coef(fit, s = fit$lambda)[-1, which.min(mse)]
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



#' Computer the consistency lower bound
#' 
#' @param df Table of (de)activation rates.
#' @param B Number of complementary pairs to draw for stability selection.
#' 
#' @import data.table
#' @import foreach

lb_fn <- function(df, B) {
  # Loop through thresholds
  lies <- function(tau) {
    # Internal consistency
    df[, dji := ifelse(drji >= tau, 1, 0)]
    df[, aji := ifelse(arji >= tau, 1, 0)]
    df[, dij := ifelse(drij >= tau, 1, 0)]
    df[, aij := ifelse(arij >= tau, 1, 0)]
    df[, int_err := ifelse((dji + aji > 1) | (dij + aij > 1), 1, 0)]
    int_err <- sum(df$int_err)
    # External consistency
    sum_ji <- df[, sum(dji + aji)]
    sum_ij <- df[, sum(dij + aij)]
    ext_err <- ifelse(min(c(sum_ji, sum_ij)) > 0, 1, 0)
    # Export
    out <- data.table('tau' = tau, 'int_err' = int_err, 'ext_err' = ext_err)
  }
  lie_df <- foreach(tt = seq_len(2 * B) / (2 * B), .combine = rbind) %do% 
    lies(tt)
  # Compute minimal thresholds
  min_int <- lie_df[int_err == 0, min(tau)]
  min_ext <- lie_df[ext_err == 0, min(tau)]
  min_two <- lie_df[int_err == 0 & ext_err == 0, min(tau)] # It's always ext
  # Export
  return(min_two)
}



#' Infer causal direction using stability selection
#' 
#' @param df Table of (de)activation rates.
#' @param lb Consistency lower bound, as computed by \code{lb_fn}.
#' @param order Causal order of interest, either \code{"ij"} or \code{"ji"}.
#' @param rule Inference rule, either \code{"R1"} or \code{"R2"}.
#' @param B Number of complementary pairs to draw for stability selection.
#' 
#' @import data.table
#' @import dplyr

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
    filter(tau > lb) %>%
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
