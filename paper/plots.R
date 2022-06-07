# Load libraries
library(data.table)
library(scales)
library(ggsci)
library(tidyverse)
library(doMC)
registerDoMC(8)


################################################################################

### Bivariate ###

# Import bivariate results
df <- readRDS('./res/lin_biv_benchmark.rds')

# Polish
df[, na := sum(g_hat == 'na'), by = .(n, g, method)]
df[, ci := sum(g_hat == 'ci'), by = .(n, g, method)]
df[, xy := sum(g_hat == 'xy'), by = .(n, g, method)]
df[, yx := sum(g_hat == 'yx'), by = .(n, g, method)]
df <- df %>%
  select(-g_hat) %>%
  pivot_longer(cols = na:yx, names_to = 'g_hat', values_to = 'prop') %>%
  select(n, g, method, g_hat, prop) %>%
  unique(.) %>%
  mutate(n = factor(paste('n =', n), 
                    #levels = c('n = 2500', 'n = 5000', 'n = 10000')),
                    levels = c('n = 5000', 'n = 10000', 'n = 20000')),
         g_hat = factor(g_hat, levels = c('xy', 'yx', 'ci', 'na'))) %>%
  as.data.table(.)
df[g == 'xy', g := '(a)']
df[g == 'ci', g := '(b)']
df[g == 'na', g := '(c)']
df[method == 'cbl', method := 'CBL']
df[method == 'constr', method := 'constraint']
df[, method := factor(method, levels = c('constraint', 'score', 'CBL'))]
df[, g_hat := recode_factor(g_hat, `xy` = 'X %->% Y', `yx` = 'Y %->% X',
                            `ci` = 'X %~% Y', `na` = 'NA')]

# Plot
p <- ggplot(df, aes(method, prop, fill = g_hat)) + 
  geom_bar(stat = 'identity') + 
  scale_fill_npg(labels = parse_format()) + 
  labs(x = 'Method', y = 'Percentage', 
       #title = 'Linear SCM',
       title = 'Nonlinear SCM',
       fill = 'Estimated\nStructure') + 
  theme_bw() +  
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position = 'none') + 
  facet_grid(g ~ n) 
#ggsave('./plots/lin_biv.pdf', width = 5, height = 6)
ggsave('./plots/nl_biv.pdf', width = 5, height = 6)
  
p + theme(legend.position = 'none')

p <- ggplot(df, aes(g_hat, prop, fill = g_hat)) + 
  geom_bar(stat = 'identity') + 
  scale_fill_npg(labels = parse_format()) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_grid(g ~ n) 

################################################################################

#### Multivariate ###

# Import data
df1 <- readRDS('./res/multivar.rds')
df2 <- readRDS('./res/multivar_rfci.rds')

# Loop through, export summary stats
extract_fn <- function(df, ni, d_zi, spi, idxi) {
  tmp <- df %>% filter(n == ni, d_z == d_zi, sp == spi, idx == idxi)
  tru <- as.numeric(tmp$amat[[1]])
  keep <- !is.na(tru)
  y <- tru[keep]
  yhat_cbl <- as.numeric(tmp$amat_cbl[[1]])[keep]
  yhat_cbl[yhat_cbl == 0.5] <- NA_real_
  yhat_ges <- as.numeric(tmp$amat_ges[[1]])[keep]
  out <- data.table(
    n = ni, d_z = d_zi, sp = spi, idx = idxi, 
    y = rep(y, times = 2), y_hat = c(yhat_cbl, yhat_ges), 
    method = rep(c('cbl', 'ges'), each = length(y))
  )
  return(out)
}
df <- foreach(nn = c(500, 1000, 2000, 4000, 8000), .combine = rbind) %:%
  foreach(dd = c(50, 100), .combine = rbind) %:%
  foreach(ss = c(1/4, 3/4), .combine = rbind) %:%
  foreach(ii = 1:20, .combine = rbind) %dopar% 
  extract_fn(df1, nn, dd, ss, ii)

method_smry <- function(m) {
  tmp <- df[method == m & !is.na(y_hat)]
  out <- tmp[, acc := sum(y == y_hat) / .N, by = .(n, d_z, sp, idx)]
  out <- unique(out[, .(n, d_z, sp, idx, method, acc)])
  return(out)
}
res1 <- foreach(mm = c('cbl', 'ges'), .combine = rbind) %dopar% 
  method_smry(mm)

# Similar but different function for RFCI
extract_fn <- function(df, ni, d_zi, spi, idxi) {
  tmp <- df %>% filter(n == ni, d_z == d_zi, sp == spi, idx == idxi)
  tru <- as.numeric(tmp$amat[[1]])
  keep <- !is.na(tru)
  y <- tru[keep]
  y_hat <- as.numeric(tmp$amat_rfci[[1]])[keep]
  out <- data.table(
    n = ni, d_z = d_zi, sp = spi, idx = idxi, 
    y, y_hat, method = 'rfci'
  )
  return(out)
}
df_a <- foreach(dd = c(50, 100), .combine = rbind) %:%
  foreach(ss = c(1/4, 3/4), .combine = rbind) %:%
  foreach(ii = 1:5, .combine = rbind) %dopar% 
  extract_fn(df2, ni = 500, dd, ss, ii)
df_b <- foreach(ss = c(1/4, 3/4), .combine = rbind) %:%
  foreach(ii = 1:5, .combine = rbind) %dopar% 
  extract_fn(df2, ni = 1000, d_zi = 50, ss, ii)
df <- rbind(df_a, df_b)
res2 <- method_smry('rfci')
res <- rbind(res1, res2)

# Polish for plotting
res[, acc_mu := mean(acc), by = .(method, n, d_z, sp)]
res[, acc_sigma := sd(acc), by = .(method, n, d_z, sp)]
df <- unique(res[, .(n, d_z, sp, method, acc_mu, acc_sigma)])
df[method == 'cbl', method := 'CBL']
df[method == 'ges', method := 'GES']
df[method == 'rfci', method := 'RFCI']
df[sp == 1/4, sparsity := 'sparsity = 0.25']
df[sp == 3/4, sparsity := 'sparsity = 0.75']

# Plot
p <- ggplot(df, aes(n, acc_mu, group = method, color = method)) + 
  geom_point(aes(shape = method), size = 3) + 
  geom_path(size = 0.5) + 
  geom_errorbar(aes(ymin = acc_mu - acc_sigma, 
                    ymax = acc_mu + acc_sigma), width = 0.05) + 
  scale_shape_manual(values = c(1, 2, 0)) +
  scale_color_npg() + 
  scale_x_log10(breaks = c(500, 1000, 2000, 4000, 8000)) + 
  #ylim(0.4, 1) + 
  labs(x = 'Sample Size', y = 'Accuracy') + 
  guides(color = guide_legend(title = 'Method'),
         shape = guide_legend(title = 'Method')) + 
  theme_bw() +
  theme(panel.border = element_rect(color = 'black', fill = NA),
        legend.position = c(0.92, 0.1), 
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.box.background = element_rect(colour = 'black')) + 
  facet_grid(d_z ~ sparsity, labeller = label_bquote(italic(d[Z])==.(d_z)))


