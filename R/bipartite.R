#' Simulated data 
#'
#' Simulated dataset of \eqn{n=200} samples with 2 foreground variables and 10
#' background variables. The design follows that of Watson & Silva (2022), with
#' \eqn{Z} drawn from a multivariate Gaussian distribution with a Toeplitz
#' covariance matrix of autocorrelation \eqn{\rho = 0.25}. Expected sparsity is
#' 0.5, signal-to-noise ratio is 2, and structural equations are linear. The 
#' ground truth for foreground variables is \eqn{X \rightarrow Y}. 
#'
#' @docType data
#'
#' @usage data(bipartite)
#'
#' @format A list with two elements: \code{x} (foreground variables), and 
#'   \code{z} (background variables).
#'
#' @keywords datasets
#'
#' @references 
#' Watson, D.S. & Silva, R. (2022). Causal discovery under a confounder blanket.
#' To appear in \emph{Proceedings of the 38th Conference on Uncertainty in 
#' Artificial Intelligence}. \emph{arXiv} preprint, 2205.05715. 
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
'bipartite'