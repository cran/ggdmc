# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

r_fastdm <- function(num_values, params, precision = 3, stop_on_error = TRUE) {
    .Call('_ggdmc_r_fastdm', PACKAGE = 'ggdmc', num_values, params, precision, stop_on_error)
}

#' Calculate likelihoods
#'
#' These function calculate likelihoods. \code{likelihood_rd} implements
#' the equations in Voss, Rothermund, and Voss (2004). These equations
#' calculate diffusion decision model (Ratcliff & Mckoon, 2008). Specifically,
#' this function implements Voss, Rothermund, and Voss's (2004) equations A1
#' to A4 (page 1217) in C++.
#'
#' @param p_vector_r a parameter vector
#' @param data data model instance
#' @param min_lik minimal likelihood.
#' @return a vector
#' @references Voss, A., Rothermund, K., & Voss, J. (2004).  Interpreting the
#' parameters of the diffusion model: An empirical validation.
#' \emph{Memory & Cognition}, \bold{32(7)}, 1206-1220. \cr\cr
#' Ratcliff, R. (1978). A theory of memory retrival. \emph{Psychological
#' Review}, \bold{85}, 238-255.
#'
#' @examples
#' model <- BuildModel(
#' p.map     = list(A = "1", B = "1", t0 = "1", mean_v = "M", sd_v = "1",
#'             st0 = "1"),
#' match.map = list(M = list(s1 = 1, s2 = 2)),
#' factors   = list(S = c("s1", "s2")),
#' constants = c(st0 = 0, sd_v = 1),
#' responses = c("r1", "r2"),
#' type      = "norm")
#'
#' p.vector <- c(A = .25, B = .35,  t0 = .2, mean_v.true = 1,
#'           mean_v.false = .25)
#' dat <- simulate(model, 1e3,  ps = p.vector)
#' dmi <- BuildDMI(dat, model)
#' den <- likelihood(p.vector, dmi)
#'
#' model <- BuildModel(
#' p.map     = list(a = "1", v = "1", z = "1", d = "1", t0 = "1", sv = "1",
#'             sz = "1", st0 = "1"),
#' constants = c(st0 = 0, d = 0),
#' match.map = list(M = list(s1 = "r1", s2 = "r2")),
#' factors   = list(S = c("s1", "s2")),
#' responses = c("r1", "r2"),
#' type      = "rd")
#'
#' p.vector <- c(a = 1, v = 1, z = 0.5, sz = 0.25, sv = 0.2, t0 = .15)
#' dat <- simulate(model, 1e2, ps = p.vector)
#' dmi <- BuildDMI(dat, model)
#' den <- likelihood (p.vector, dmi)
#'
#' @export
likelihood <- function(p_vector_r, data, min_lik = 1e-10) {
    .Call('_ggdmc_likelihood', PACKAGE = 'ggdmc', p_vector_r, data, min_lik)
}

p_df <- function(p_vector_r, cell, mtype, pnames, parnames, dim0, dim1, dim2, allpar, model, isr1, n1idx, n1order) {
    .Call('_ggdmc_p_df', PACKAGE = 'ggdmc', p_vector_r, cell, mtype, pnames, parnames, dim0, dim1, dim2, allpar, model, isr1, n1idx, n1order)
}

#' Generate Random Responses of the LBA Distribution
#'
#' \code{rlba_norm} used the LBA process to generate response times and
#' responses.
#'
#' @param n is the numbers of observation.
#' @param A start point upper bound, a vector of a scalar.
#' @param b decision threshold, a vector or a scalar.
#' @param mean_v mean drift rate vector
#' @param sd_v standard deviation of drift rate vector
#' @param t0 non-decision time, a vector.
#' @param st0 non-decision time variation, a vector.
#' @param posdrift if exclude negative drift rates
#'
#' @return a n x 2 matrix of response time (first column) and responses (second
#' column).
#' @export
rlba_norm <- function(n, A, b, mean_v, sd_v, t0, st0, posdrift) {
    .Call('_ggdmc_rlba_norm', PACKAGE = 'ggdmc', n, A, b, mean_v, sd_v, t0, st0, posdrift)
}

rprior_mat <- function(prior, n) {
    .Call('_ggdmc_rprior_mat', PACKAGE = 'ggdmc', prior, n)
}

init_mcmc <- function(data_or_samples, prior, nchain, nmc, thin, report, rp, gammamult, pm_Hu, pm_BT, block, add = FALSE, is_old = FALSE) {
    .Call('_ggdmc_init_mcmc', PACKAGE = 'ggdmc', data_or_samples, prior, nchain, nmc, thin, report, rp, gammamult, pm_Hu, pm_BT, block, add, is_old)
}

init_hier <- function(data_or_samples, prior, lprior, sprior, nchain, nmc, thin, report, rp, gammamult, pm_Hu, pm_BT, block, add = FALSE, is_old = FALSE) {
    .Call('_ggdmc_init_hier', PACKAGE = 'ggdmc', data_or_samples, prior, lprior, sprior, nchain, nmc, thin, report, rp, gammamult, pm_Hu, pm_BT, block, add, is_old)
}

#' Truncated Normal Distribution
#'
#' Random number generation, probability density and cumulative density
#' functions for truncated normal distribution.
#'
#' @param x,q vector of quantiles;
#' @param n number of observations. n must be a scalar.
#' @param p1 mean (must be scalar).
#' @param p2 standard deviation (must be scalar).
#' @param lower lower truncation value (must be scalar).
#' @param upper upper truncation value (must be scalar).
#' @param lt lower tail. If TRUE (default) probabilities are \code{P[X <= x]},
#' otherwise, \code{P[X > x]}.
#' @param lg log probability. If TRUE (default is FALSE) probabilities p are
#' given as \code{log(p)}.
#' @return a column vector.
#' @examples
#' ## rtn example
#' dat1 <- rtnorm(1e5, 0, 1, 0, Inf)
#' hist(dat1, breaks = "fd", freq = FALSE, xlab = "",
#'      main = "Truncated normal distributions")
#'
#' ## dtn example
#' x <- seq(-5, 5, length.out = 1e3)
#' dat1 <- dtnorm(x, 0, 1, -2, 2, 0)
#' plot(x, dat1, type = "l", lwd = 2, xlab = "", ylab= "Density",
#'      main = "Truncated normal distributions")
#'
#' ## ptn example
#' x <- seq(-10, 10, length.out = 1e2)
#' mean <- 0
#' sd <- 1
#' lower <- 0
#' upper <- 5
#' dat1 <- ptnorm(x, 0, 1, 0, 5, lg = TRUE)
#' @export
dtnorm <- function(x, p1, p2, lower, upper, lg = FALSE) {
    .Call('_ggdmc_dtnorm', PACKAGE = 'ggdmc', x, p1, p2, lower, upper, lg)
}

#' @rdname dtnorm
#' @export
rtnorm <- function(n, p1, p2, lower, upper) {
    .Call('_ggdmc_rtnorm', PACKAGE = 'ggdmc', n, p1, p2, lower, upper)
}

#' @rdname dtnorm
#' @export
ptnorm <- function(q, p1, p2, lower, upper, lt = TRUE, lg = FALSE) {
    .Call('_ggdmc_ptnorm', PACKAGE = 'ggdmc', q, p1, p2, lower, upper, lt, lg)
}

