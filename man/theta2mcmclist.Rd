% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analysis.R
\name{theta2mcmclist}
\alias{theta2mcmclist}
\alias{phi2mcmclist}
\title{Convert theta to a mcmc List}
\usage{
theta2mcmclist(
  x,
  start = 1,
  end = NA,
  split = FALSE,
  subchain = FALSE,
  nsubchain = 3,
  thin = NA
)

phi2mcmclist(
  x,
  start = 1,
  end = NA,
  split = FALSE,
  subchain = FALSE,
  nsubchain = 3
)
}
\arguments{
\item{x}{posterior samples}

\item{start}{start iteration}

\item{end}{end iteraton}

\item{split}{whether to divide one MCMC sequence into two sequences.}

\item{subchain}{boolean swith convert only a subset of chains}

\item{nsubchain}{indicate the number of chains in the subset}

\item{thin}{thinning lenght of the posterior samples}
}
\description{
Extracts the parameter array (ie theta) from posterior samples of a
partiipant and convert it to a \pkg{coda} mcmc.list.
}
\details{
\code{phi2mcmclist} extracts the phi parameter array, which stores
the location and scale parameters at the hyper level.
}
\examples{
\dontrun{
model <- BuildModel(
p.map     = list(a = "RACE", v = c("S", "RACE"), z = "RACE", d = "1",
            sz = "1", sv = "1", t0 = c("S", "RACE"), st0 = "1"),
match.map = list(M = list(gun = "shoot", non = "not")),
factors   = list(S = c("gun", "non"), RACE = c("black", "white")),
constants = c(st0 = 0, d = 0, sz = 0, sv = 0),
responses = c("shoot", "not"),
type      = "rd")

pnames <- GetPNames(model)
npar   <- length(pnames)
pop.mean  <- c(1, 1, 2.5, 2.5, 2.5, 2.5, .50, .50, .4, .4, .4, .4)
pop.scale <- c(.15, .15, 1, 1, 1, 1, .05, .05, .05, .05, .05, .05)
names(pop.mean)  <- pnames
names(pop.scale) <- pnames
pop.prior <- BuildPrior(
  dists = rep("tnorm", npar),
  p1    = pop.mean,
  p2    = pop.scale,
  lower = c(rep(0, 2), rep(-5, 4), rep(0, 6)),
  upper = c(rep(5, 2), rep(7, 4), rep(2, 6)))
p.prior <- BuildPrior(
  dists = rep("tnorm", npar),
  p1    = pop.mean,
  p2    = pop.scale*10,
  lower = c(rep(0, 2), rep(-5, 4), rep(0, 6)),
  upper = c(rep(10, 2), rep(NA, 4), rep(5, 6)))
mu.prior <- BuildPrior(
  dists = rep("tnorm", npar),
  p1    = pop.mean,
  p2    = pop.scale*10,
  lower = c(rep(0,  2), rep(-5, 4), rep(0, 6)),
  upper = c(rep(10, 2), rep(NA, 4), rep(5, 6)))
sigma.prior <- BuildPrior(
  dists = rep("beta", npar),
  p1    = rep(1, npar),
  p2    = rep(1, npar),
  upper = rep(2, npar))
names(sigma.prior) <- GetPNames(model)
priors <- list(pprior=p.prior, location=mu.prior, scale=sigma.prior)
dat    <- simulate(model, nsim = 10, nsub = 10, prior = pop.prior)
dmi    <- BuildDMI(dat, model)
ps     <- attr(dat, "parameters")

fit0 <- StartNewsamples(dmi, priors)
fit  <- run(fit0)

tmp1 <- theta2mcmclist(fit[[1]])
tmp2 <- theta2mcmclist(fit[[2]], start = 10, end = 90)
tmp3 <- theta2mcmclist(fit[[3]], split = TRUE)
tmp4 <- theta2mcmclist(fit[[4]], subchain = TRUE)
tmp5 <- theta2mcmclist(fit[[5]], subchain = TRUE, nsubchain = 4)
tmp6 <- theta2mcmclist(fit[[6]], thin = 2)
}

}
