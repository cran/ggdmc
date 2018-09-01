#### Automatic Sampling ---------------------------------------
##' Four pMCMC checking tools for automatic sampling
##'
##' The function tests whether MCMC chains encounter a parameter region
##' difficult to search (ie get stuck):
##'
##' @param x posterior samples
##' @param cut the criteria for suggesting abnormal chains found
##' @param verbose print more information
##'
##' @export
isstuck <- function(x, cut = 10, verbose = FALSE) {
  if (verbose) cat("Stuck chains check\n")
  stucks <- PickStuck(x, cut = cut, verbose = verbose)
  fail <- length( stucks != 0)
  if (verbose) {
    if (!fail) cat(": OK\n") else cat(paste(":", length(stucks),"\n"))
  }
  return(fail)
}


gmcmc <- function(x) {
  coda::mcmc(matrix(aperm(x$theta, c(1, 3, 2)),
                    ncol = dim(x$theta)[2], dimnames = list(NULL, dimnames(x$theta)[[2]])))
}



##' Four pMCMC checking tools for automatic sampling
##'
##' The function tests whether MCMC chains converge prematurelly:
##'
##' @param x posterior samples
##' @param p1 the range of the head of MCMC chains
##' @param p2 the range of the tail of the MCMC chains
##' @param cut_location how far away a location chains been considered as stuck
##' @param cut_scale how far away a scale chains been considered as stuck
##' @param verbose print more information
##' @export
##' @importFrom stats median
##' @importFrom stats IQR
##' @export
isflat <- function(x, p1 = 1/3, p2 = 1/3, cut_location = 0.25,
  cut_scale = Inf, verbose = FALSE) {

  mat  <- gmcmc(x)
  xlen <- round(dim(mat)[1] * p1)
  ylen <- round(dim(mat)[1] * p2)
  ## change in mean relative to robust SD
  m.zs <- apply(mat, 2, function(xx) {
    m1 <- median(xx[1:xlen])
    m2 <- median(xx[(length(xx)-ylen):length(xx)])
    abs(m1 - m2) / IQR(xx)
  })

  names(m.zs) <- paste("m", names(m.zs), sep="_")
  fail <- any(m.zs > cut_location)
  if (!fail) out <- "" else
    out <- paste(names(m.zs)[m.zs==max(m.zs)], "=", round(max(m.zs), 2))

  if ( is.finite(cut_scale) ) {
    ## Change ini IQR reltiave to overall IQR
    s.zs <- apply(mat, 2, function(xx){
      m1 <- IQR(xx[1:xlen])
      m2 <- IQR(xx[(length(xx)-ylen):length(xx)])
      abs(m1 - m2)/IQR(xx)
    })
    names(s.zs) <- paste("s",names(s.zs), sep = "_")
    if (out != "") out <- paste(out,", ", sep = "")
    if (any(s.zs > cut_scale)) out <-
      paste(out, names(s.zs)[s.zs == max(s.zs)], "=", round(max(s.zs),2))
    fail <- fail | any(s.zs > cut_scale)
  }

  if (verbose) {
    cat("Flat check\n")
    print(round(m.zs, 2))
    if ( is.finite(cut_scale) )
      print(round(s.zs,2)) else
        if (!fail) cat(": OK\n") else
          cat(paste(":",out,"\n"))
  }
  fail
}

##' Four pMCMC checking tools for automatic sampling
##'
##' The function tests whether MCMC chains mixed well.
##'
##' @param x posterior samples
##' @param cut psrf criterion for well mixed
##' @param split whether to split MCMC chains. This is an argument passing to
##' gelman function
##' @param verbose print more information
##' @seealso \code{\link{gelman}})
##' @export
ismixed <- function(x, cut = 1.01, split = TRUE, verbose = FALSE) {

  tmp <- gelman(x, split = split)
  gds <- c(tmp$mpsrf, tmp$psrf[,1])
  fail <- max(gds) > cut
  if (verbose) {
    cat ("Mixing check\n")
    print(round(gds, 2))
    if (!fail) {
      cat(": OK\n")
    } else {
      nam <- names(gds)[gds==max(gds)]
      cat (paste(":", nam, "=", round(max(gds), 2), "\n"))
    }
  }
  fail
}

##' Four pMCMC checking tools for automatic sampling
##'
##' The function tests whether MCMC chains have drawn enough samples.
##'
##' @param x posterior samples
##' @param minN minimal effective sample sizes
##' @param nfun mean or median function
##' @param verbose print more information
##' @export
iseffective <- function(x, minN, nfun, verbose = FALSE) {
  n <- do.call(nfun, list(effectiveSize(x, verbose = verbose)))
  fail <- n < minN
  if (verbose) {
    cat("Length check")
    if (!fail) cat(": OK\n") else cat(paste(":",n,"\n"))
  }
  fail
}

### MCMC -------------------------------------------------------
##' Convert theta to a mcmc List
##'
##' Extracts the parameter array (ie theta) from posterior samples of a
##' partiipant and convert it to a \pkg{coda} mcmc.list.
##'
##' \code{phi2mcmclist} extracts the phi parameter array, which store
##' the location and scale parameters at the hyper level.
##'
##' @param x posterior samples
##' @param start start iteration
##' @param end end iteraton
##' @param split whether to divide one MCMC sequence into two sequences.
##' @param subchain boolean swith convert only a subset of chains
##' @param nsubchain indicate the number of chains in the subset
##' @param thin thinning lenght of the posterior samples
##' @importFrom coda mcmc mcmc.list
##' @examples
##' \dontrun{
##' model <- BuildModel(
##' p.map     = list(a = "RACE", v = c("S", "RACE"), z = "RACE", d = "1", sz = "1",
##'   sv = "1", t0 = c("S", "RACE"), st0 = "1"),
##' match.map = list(M = list(gun = "shoot", non = "not")),
##' factors   = list(S = c("gun", "non"), RACE = c("black", "white")),
##' constants = c(st0 = 0, d = 0, sz = 0, sv = 0),
##' responses = c("shoot", "not"),
##' type      = "rd")
##'
##' pnames <- GetPNames(model)
##' npar <- length(pnames)
##' pop.mean  <- c(1, 1, 2.5, 2.5, 2.5, 2.5, .50, .50, .4, .4, .4, .4)
##' pop.scale <- c(.15, .15, 1, 1, 1, 1, .05, .05, .05, .05, .05, .05)
##' names(pop.mean)  <- pnames
##' names(pop.scale) <- pnames
##' pop.prior <- BuildPrior(
##'   dists = rep("tnorm", npar),
##'   p1    = pop.mean,
##'   p2    = pop.scale,
##'   lower = c(rep(0, 2), rep(-5, 4), rep(0, 6)),
##'   upper = c(rep(5, 2), rep(7, 4), rep(2, 6)))
##' p.prior <- BuildPrior(
##'   dists = rep("tnorm", npar),
##'   p1    = pop.mean,
##'   p2    = pop.scale*10,
##'   lower = c(rep(0, 2), rep(-5, 4), rep(0, 6)),
##'   upper = c(rep(10, 2), rep(NA, 4), rep(5, 6)))
##' mu.prior <- BuildPrior(
##'   dists = rep("tnorm", npar),
##'   p1    = pop.mean,
##'   p2    = pop.scale*10,
##'   lower = c(rep(0,  2), rep(-5, 4), rep(0, 6)),
##'   upper = c(rep(10, 2), rep(NA, 4), rep(5, 6)))
##' sigma.prior <- BuildPrior(
##'   dists = rep("beta", npar),
##'   p1    = rep(1, npar),
##'   p2    = rep(1, npar),
##'   upper = rep(2, npar))
##' names(sigma.prior) <- GetPNames(model)
##' pp.prior <- list(mu.prior, sigma.prior)
##'
##' dat <- simulate(model, nsim = 30, nsub = 10, p.prior = pop.prior)
##' dmi <- BuildDMI(dat, model)
##' ps <- attr(dat, "parameters")
##'
##' hsam <- run(StartNewHypersamples(1e2, dmi, p.prior, pp.prior, 1),
##'   pm = .1, hpm = .1, report = 20)
##'
##' tmp1 <- theta2mcmclist(hsam[[1]])
##' tmp2 <- theta2mcmclist(hsam[[2]], start = 10, end = 90)
##' tmp3 <- theta2mcmclist(hsam[[3]], split = TRUE)
##' tmp4 <- theta2mcmclist(hsam[[4]], subchain = TRUE)
##' tmp5 <- theta2mcmclist(hsam[[5]], subchain = TRUE, nsubchain = 4)
##' tmp6 <- theta2mcmclist(hsam[[6]], thin = 2)
##' }
##'
##' @export
theta2mcmclist <- function(x, start = 1, end = NA, split = FALSE,
  subchain = FALSE, nsubchain = 3, thin = NA) {

  if (is.na(thin)) thin <- x$thin
  if (is.na(end)) end <- x$nmc
  nchain <- x$n.chains

  if (subchain) {
    message("pMCMC diagnosis randomly select a subset of chains: ", appendLF = FALSE)
    cidx <- base::sample(1:nchain, nsubchain)
    cat(cidx, "\n")
    nchain <- nsubchain
  } else {
    cidx <- 1:nchain
  }

  lst <- vector(mode = "list", length = nchain * ifelse(split, 2, 1))
  iter <- start:end

  if (split) is.in <- !as.logical(iter %% 2) else is.in <- rep(TRUE, length(iter))

  if (split) {
    not.is.in <- !is.in
    if ( sum(not.is.in) > sum(is.in) ) {
      not.is.in[1:2][not.is.in[1:2]] <- FALSE
    } else if ( sum(not.is.in) < sum(is.in) ) {
      is.in[1:2][is.in[1:2]] <- FALSE
    }
  }

  for (i in 1:nchain) {
    lst[[i]] <- coda::mcmc( t(x$theta[cidx[i], , iter[is.in]]),
      thin = thin)
  }

  if (split) {
    for (i in 1:nchain) {
      lst[[i + nchain]] <- coda::mcmc( t(x$theta[cidx[i], , iter[not.is.in]]),
        thin = thin)
    }
  }

  attr(lst, "nchain") <- nchain
  attr(lst, "npar")   <- x$n.pars
  attr(lst, "thin")   <- thin
  attr(lst, "iter")   <- iter
  attr(lst, "pnames") <- x$p.names
  attr(lst, "nmc")    <- x$nmc
  attr(lst, "start")  <- start
  attr(lst, "end")    <- end
  return(lst)
}

##' @rdname theta2mcmclist
##' @importFrom coda mcmc mcmc.list
##' @export
phi2mcmclist <- function(x, start = 1, end = NA, split = FALSE,
  subchain = FALSE, nsubchain = 3) {

  thin   <- x$thin   ## x == hyper
  nchain <- x$n.chains

  if (subchain) {
    message("pMCMC diagnosis randomly select a subset of chains: ", appendLF = FALSE)
    chain.idx <- base::sample(1:nchain, nsubchain)
    cat(chain.idx, "\n")
    nchain <- nsubchain
  } else {
    chain.idx <- 1:nchain
  }

  # Parameters that are not constants
  ok1 <- lapply(x$pp.prior,function(x){
    lapply(x,function(y){attr(y, "dist") != "constant"})})
  ok2 <- paste(names(ok1[[2]])[unlist(ok1[[2]])], "h2", sep=".")
  ok1 <- paste(names(ok1[[1]])[unlist(ok1[[1]])], "h1", sep=".")
  if ( is.na(end) ) end <- x$nmc

  lst <- vector(mode = "list", length = nchain)
  indx <- start:end

  if (split) is.in <- !as.logical(indx %% 2) else is.in <- rep(TRUE, length(indx))
  if (split) {
    not.is.in <- !is.in
    if ( sum(not.is.in) > sum(is.in) ) {
      not.is.in[1:2][not.is.in[1:2]] <- FALSE
    } else if ( sum(not.is.in) < sum(is.in) ) {
      is.in[1:2][is.in[1:2]] <- FALSE
    }
  }

  for (i in 1:nchain) {
    tmp1 <- t(x$phi[[1]][chain.idx[i], , indx[is.in]]) ## nmc x npar matrix
    ## attach parnames with h1
    dimnames(tmp1)[[2]] <- paste(dimnames(tmp1)[[2]],"h1", sep=".")
    tmp1 <- tmp1[,ok1]  ## exclude constant parameter

    ## same thing on scale
    tmp2 <- t(x$phi[[2]][chain.idx[i], , indx[is.in]])
    dimnames(tmp2)[[2]] <- paste(dimnames(tmp2)[[2]],"h2", sep=".")
    tmp2 <- tmp2[,ok2]

    ## Remove cases with !has.sigma
    tmp2 <- tmp2[,!apply(tmp2, 2, function(x){all(is.na(x))})]
    lst[[i]] <- coda::mcmc(cbind(tmp1, tmp2), thin = thin)
  }

  if (split) {
    for (i in 1:nchain) {
      tmp1 <- t(x$phi[[1]][chain.idx[i], , indx[not.is.in]])
      dimnames(tmp1)[[2]] <- paste(dimnames(tmp1)[[2]],"h1", sep=".")
      tmp1 <- tmp1[,ok1]
      tmp2 <- t(x$phi[[2]][chain.idx[i], , indx[not.is.in]])
      dimnames(tmp2)[[2]] <- paste(dimnames(tmp2)[[2]],"h2", sep=".")
      tmp2 <- tmp2[,ok2]
      # Remove cases with !has.sigma
      tmp2 <- tmp2[,!apply(tmp2,2,function(x){all(is.na(x))})]
      lst[[i + nchain]] <- coda::mcmc(cbind(tmp1,tmp2), thin = thin)
    }
  }

  attr(lst, "nchain") <- nchain
  attr(lst, "npar")   <- x$n.pars * 2
  attr(lst, "thin")   <- thin
  attr(lst, "iter")   <- indx
  attr(lst, "pnames") <- dimnames(lst[[1]])[[2]]
  attr(lst, "nmc")    <- x$nmc
  attr(lst, "start")  <- start
  attr(lst, "end")    <- end
  return(lst)
}


##' Potential scale reduction factor
##'
##' \code{gelman} calls \pkg{coda} gelman.diag to get R hats for one
##' or a list of subjects. It calculates at the either data or hyper level.
##'
##' @param x posterior samples
##' @param hyper a Boolean switch, indicating posterior samples are from
##' hierarchical modeling
##' @param start start iteration
##' @param end end iteration
##' @param confidence confident inteval
##' @param transform turn on transform
##' @param autoburnin turn on auto burnin
##' @param multivariate multivariate Boolean switch
##' @param split split whether split mcmc chains; When split is TRUE, the function
##' doubles the number of chains by spliting into 1st and 2nd halves.
##' @param subchain whether only calculate a subset of chains
##' @param nsubchain indicate how many chains in a subset
##' @param verbose print more information
##' @param digits print out how many digits
##' @param ... arguments passing to \code{coda} gelman.diag.
##' @importFrom coda gelman.diag
##' @export
##' @examples
##' \dontrun{
##' rhat1 <- hgelman(hsam); rhat1
##' rhat2 <- hgelman(hsam, end = 51); rhat2
##' rhat3 <- hgelman(hsam, confidence = .90); rhat3
##' rhat4 <- hgelman(hsam, transform = FALSE); rhat4
##' rhat5 <- hgelman(hsam, autoburnin = TRUE); rhat5
##' rhat6 <- hgelman(hsam, split = FALSE); rhat6
##' rhat7 <- hgelman(hsam, subchain = TRUE); rhat7
##' rhat8 <- hgelman(hsam, subchain = TRUE, nsubchain = 4);
##' rhat9 <- hgelman(hsam, subchain = TRUE, nsubchain = 4,
##' digits = 1, verbose = TRUE);
##'
##' hat1 <- gelman(hsam[[1]], multivariate = FALSE); hat1
##' hat2 <- gelman(hsam[[1]], hyper = TRUE, verbose = TRUE); hat2
##' hat3 <- gelman(hsam, hyper = TRUE, verbose = TRUE); hat3
##' hat4 <- gelman(hsam, multivariate = TRUE, verbose = FALSE);
##' hat5 <- gelman(hsam, multivariate = FALSE, verbose = FALSE);
##' hat6 <- gelman(hsam, multivariate = FALSE, verbose = TRUE);
##' hat7 <- gelman(hsam, multivariate = T, verbose = TRUE);
##' }
gelman <- function(x, hyper = FALSE, start = 1, end = NA, confidence = 0.95,
  transform=TRUE, autoburnin = FALSE, multivariate = TRUE, split = TRUE,
  subchain = FALSE, nsubchain = 3, digits = 2, verbose = FALSE, ...) {

  if ( hyper ) {
    if (verbose) message("Diagnosing the hyper parameters, phi")
    hyper <- attr(x, "hyper")
    if (is.null(hyper)) stop("Posterior samples are not from a hierarchical fit")
    thin <- hyper$thin
    if (is.na(end)) end <- hyper$nmc

    mcmclist <- phi2mcmclist(hyper, start, end, split, subchain, nsubchain)
    out <- coda::gelman.diag(mcmclist, confidence, transform, autoburnin,
      multivariate)

  } else {
    ## if x is one subject samples, we should found an elemnet called theta
    if ( !is.null(x$theta) ) {
      if (verbose) message("Diagnosing a single participant, theta")
      if (is.na(end)) end <- x$nmc
      mcmclist <- theta2mcmclist(x, start, end, split, subchain, nsubchain)
      out <- coda::gelman.diag(mcmclist, confidence, transform, autoburnin,
        multivariate)
    } else {
      if (verbose) message("Diagnosing theta for many participants separately")
      out <- lapply(x, function(xx) {
        step1 <- theta2mcmclist(xx, start, end, split, subchain, nsubchain)
        step2 <- coda::gelman.diag(step1, confidence, transform, autoburnin,
          multivariate)
        return(step2)
      })

      names(out)  <- names(x)
      tmp <- unlist(lapply(out, function(x){x$mpsrf}))

      if (verbose) {
        if (!multivariate) stop("Must set multivariate to TRUE")
        printthis <- c(mean(tmp), sort(tmp))
        names(printthis) <- c("Mean", names(sort(tmp)))
        print( round(printthis, digits) )
      }
    }
  }

  return(out)
}


##' @rdname gelman
##' @export
hgelman <- function(x, start = 1, end = NA, confidence = 0.95, transform = TRUE,
  autoburnin = FALSE, split = TRUE, subchain = FALSE,
  nsubchain = 3, digits = 2, verbose = FALSE, ...) {

  step1 <- lapply(gelman(x, start = start, end = end, confidence = confidence,
    transform = transform, autoburnin = autoburnin, multivariate = TRUE,
    split=split, subchain = subchain, nsubchain = nsubchain, verbose = FALSE),
    function(x) {x$mpsrf} )

  snames <- names(sort(unlist(step1)))
  out <- sort(unlist(step1)) ## non-hyper

  if ( any(names(attributes(x)) == "hyper") ) {
    hyper_gd <- gelman(x, hyper = TRUE, start = start, end = end,
      confidence = confidence, transform = transform, autoburnin = autoburnin,
      multivariate = TRUE, split = split, subchain = subchain,
      nsubchain = nsubchain, verbose = FALSE)

    out <- c(hyper_gd$mpsrf, out)
    names(out) <- c("hyper", snames)
  }

  if (verbose) print(round(out, digits))
  invisible(out)
}

gelman_mpsrf <- function(mcmclist, autoburnin, transform) {

  ## robust version ONLY USED IN sampling.R IN run.converge.dmc
  ## SHOULD BE ROLLED OUT OVER FOLLOWING FUNCITONS TO AVOID CRASHES OF
  ## AUTO PROCEDURES.

  gd <- try(gelman.diag(mcmclist,
    autoburnin = autoburnin, transform = transform), silent = TRUE)
  if (class(gd)=="try-error") Inf else gd$mpsrf
}

##' @rdname effectiveSize
##' @importFrom coda effectiveSize
##' @export
effectiveSize_hyper <- function(x, start, end, digits, verbose) {
  hyper <- attr(x, "hyper")
  if (is.na(end)) end <- hyper$nmc
  phimcmc <- phi2mcmclist(hyper, start = start, end = end)
  out <- coda::effectiveSize(phimcmc)
  if (verbose) print(round(out, digits))
  invisible(return(out))
}

##' @rdname effectiveSize
##' @importFrom coda effectiveSize
##' @importFrom stats sd
##' @export
effectiveSize_many <- function(x, start, end, verbose) {
  out <- lapply(x, function(xx) {
    if (is.na(end)) end <- xx$nmc
    coda::effectiveSize(theta2mcmclist(xx, start, end))
  })

  if (verbose) {
    p1 <- round(apply(data.frame(out), 1, mean))
    p2 <- round(apply(data.frame(out), 1, sd))
    p3 <- round(apply(data.frame(out), 1, max))
    p4 <- round(apply(data.frame(out), 1, min))
    print_out <- rbind(p1, p2, p3, p4)
    rownames(print_out) <- c("MEAN", "SD", "MAX", "MIN")
    print(print_out)
  }
  invisible(return(out))
}

##' @rdname effectiveSize
##' @importFrom coda effectiveSize
##' @export
effectiveSize_one <- function(x, start, end, digits, verbose) {
  if (is.na(end)) end <- x$nmc
  out <- coda::effectiveSize(theta2mcmclist(x, start = start, end = end))
  if (verbose) print(round(out, digits))
  invisible(return(out))
}

##' Effective sample size
##'
##' \code{effectiveSize} calls \pkg{coda} effectiveSize to calculate
##' effective posterior sample size.
##'
##' @param x a samples object
##' @param hyper a switch to extract hyper attribute and calculate it
##' @param start starting iteration
##' @param end ending iteraton
##' @param digits printing digits
##' @param verbose printing more information
##' @export
##' @examples
##' #################################40
##' ## effectiveSize example
##' #################################40
##' \dontrun{
##' es1 <- effectiveSize_one(hsam[[1]], 1, 100, 2, TRUE)
##' es2 <- effectiveSize_one(hsam[[1]], 1, 100, 2, FALSE)
##' es3 <- effectiveSize_many(hsam, 1, 100, TRUE)
##' es4 <- effectiveSize_many(hsam, 1, 100, FALSE)
##' es5 <- effectiveSize_hyper(hsam, 1, 100, 2)
##' es6 <- effectiveSize(hsam, TRUE, 1, 100, 2, TRUE)
##' es7 <- effectiveSize(hsam, TRUE, 1, 100, 2, FALSE)
##' es8 <- effectiveSize(hsam, FALSE, 1, 100, 2, TRUE)
##' es9 <- effectiveSize(hsam, FALSE, 1, 100, 2, FALSE)
##' es10 <- effectiveSize(hsam[[1]], FALSE, 1, 100, 2, TRUE)
##' }
effectiveSize <- function(x, hyper = FALSE, start = 1, end = NA,
  digits = 0, verbose = FALSE) {
  if (hyper) {
    out <- effectiveSize_hyper(x, start, end, digits, verbose)
  } else if (!is.null(x$theta)){
    out <- effectiveSize_one(x, start, end, digits, verbose)
  } else {
    out <- effectiveSize_many(x, start, end, verbose)
  }
}


### Summary ------------------------------------------------------
##' @importFrom coda spectrum0.ar
safespec0 <- function(x) {
  result <- try(coda::spectrum0.ar(x)$spec)
  if (class(result) == "try-error") result <- NA
  if (class(result) == "try") result <- NA
  result
}

##' Summary statistic for posterior samples
##'
##' Calculate summary statistics for pMCMC posterior samples
##'
##' @param object posterior samples
##' @param prob summary quantile summary
##' @param ... other arguments passing in
##' @export
summary_mcmc_list <- function(object, prob = c(0.025, 0.25, 0.5, 0.75, 0.975),
  ...) {

  nchain <- attr(object, "nchain")
  npar   <- attr(object, "npar")
  thin   <- attr(object, "thin")
  pnames <- attr(object, "pnames")
  nmc    <- attr(object, "nmc")
  start  <- attr(object, "start")
  end    <- attr(object, "end")
  iter   <- attr(object, "iter")
  niter  <- length(iter)

  statnames <- c("Mean", "SD", "Naive SE", "Time-series SE")
  varstats <- matrix(nrow = npar, ncol = length(statnames),
    dimnames = list(pnames, statnames))
  xtsvar <- matrix(nrow = nchain, ncol = npar)

  if (is.matrix(object[[1]])) {
    for (i in 1:nchain) {
      for (j in 1:npar) {
        xtsvar[i, j] <- safespec0(object[[i]][, j])
      }
      xlong <- do.call("rbind", object)
    }
  } else {
    for (i in 1:nchain) xtsvar[i, ] <- safespec0(object[[i]])
    xlong <- as.matrix(object)
  }

  xmean    <- matrixStats::colMeans2(xlong)
  xvar     <- matrixStats::colVars(xlong)
  xtsvar2  <- matrixStats::colMeans2(xtsvar)
  varquant <- matrixStats::colQuantiles(xlong, probs = prob)

  varstats[, 1] <- xmean
  varstats[, 2] <- sqrt(xvar)
  varstats[, 3] <- sqrt(xvar/niter * nchain)
  varstats[, 4] <- sqrt(xtsvar2/niter * nchain)
  varquant <- drop(varquant)
  varstats <- drop(varstats)
  out <- list(statistics = varstats, quantiles = varquant,
    start = start, end = end, thin = thin, nchain = nchain)
  return(out)
}

summary_hyper <- function(x, start, end, hmeans, hci, prob, digits, verbose) {

  if (verbose) message("Summarise hierarchical model")
  hyper <- attr(x, "hyper")
  if (is.null(hyper)) stop("Samples are not from a hierarhcial model fit")
  # message("end is missing detected.")
  if (is.na(end)) end <- hyper$nmc
  npar <- hyper$n.pars
  hest <- summary_mcmc_list(phi2mcmclist(hyper, start, end), prob = prob)

  if (hmeans) {
    h1 <- hest$statistics[1:npar, "Mean"]
    h2 <- hest$statistics[(1 + npar):(2 * npar), "Mean"]
    out <- round(rbind(h1, h2), digits)
    colnames(out) <- hyper$p.names

  } else if (hci) {
    quan <- hest$quantiles[, prob]
    conf <- cbind( quan[1:npar, ], quan[(1 + npar):(2 * npar), ])
    parname_noh <- unlist(strsplit(dimnames(conf)[[1]], ".h1"))
    rep_percent <- dimnames(conf)[[2]]
    per_names <- paste(rep(c("L", "S"), each = length(prob)), colnames(conf))
    dimnames(conf) <- list(parname_noh, per_names)
    out <- round(conf, digits)
  } else {
    out <- hest
  }

  return(out)
}

summary_one <- function(x, start, end, prob, verbose) {
  if (verbose) message("Single Participant")
  if (is.na(end)) end <- x$nmc
  return(summary_mcmc_list(theta2mcmclist(x, start, end),
    prob = prob))
}

##' @importFrom matrixStats colMeans2
summary_many <- function(x, start, end, prob, verbose) {
  if(verbose) message("Summary each participant separately")
  out1 <- lapply(x, function(xx, starti, endi, probi) {
          step1 <- summary_mcmc_list(theta2mcmclist(xx, starti, endi),
            prob = probi)
          return(step1)
      }, start, end, prob)
  names(out1) <- names(x)

  if (verbose) {
    return(out1)
  } else {
    df_form <- t(data.frame(lapply(out1, function(xx){xx[[1]][, 1]})))
    out2 <- rbind(df_form, matrixStats::colMeans2(df_form))
    row.names(out2) <- c(names(x), "Mean")
    return(out2)
  }
}

summary_recoverone <- function(object, start, end, ps, digits, prob, verbose) {

  if (missing(start)) start <- object$start
  if (missing(end)) end <- object$end

  qs <- summary_one(object, start, end, prob, FALSE)$quantiles
  parnames <- dimnames(qs)[[1]]

  if (!is.null(ps) && (!all(parnames %in% names(ps))))
    stop("Names of p.vector do not match parameter names in samples")

  est  <- qs[names(ps), "50%"]
  op.vector <- ps[order(names(ps))]
  oest <- est[order(names(est))]
  bias <- oest- op.vector

  lo  <- qs[names(ps), "2.5%"]
  hi  <- qs[names(ps), "97.5%"]
  olo <- lo[order(names(lo))]
  ohi <- hi[order(names(hi))]

  out  <- rbind(
    'True'          = op.vector,
    '2.5% Estimate' = olo,
    '50% Estimate'  = oest,
    '97.5% Estimate'= ohi,
    'Median-True'   = bias)

  if (verbose) print(round(out, digits))
  invisible(return(out))
}

summary_recovermany <- function(object, start, end, ps, digits, prob) {

  est <- summary_many(object, start, end, prob, TRUE)

  df_form <- t(data.frame(lapply(est, function(x){x[[1]][, 1]})))

  mean.est <- matrixStats::colMeans2(df_form)
  mean.ps <- matrixStats::colMeans2(ps)
  sd.est <- matrixStats::colSds(df_form)
  sd.ps <- matrixStats::colSds(ps)

  pnames <- colnames(ps)
  loc <- rbind(mean.est, mean.ps, mean.ps - mean.est)
  sca <- rbind(sd.est, sd.ps, sd.ps - sd.est)
  out <- rbind(loc, sca)

  rownames(out) <- c("Mean", "True", "Diff", "Sd", "True", "Diff")
  colnames(out) <- object[[1]]$p.names
  print(round(out, digits))
  invisible(return(out))
}

summary_recoverhyper <- function(object, start, end, ps, type, digits, prob,
  verbose) {

  hyper <- attr(object, "hyper")
  samples <- list(theta = hyper$phi[[type]])
  samples$n.chains <- hyper$n.chains
  samples$nmc  <- hyper$nmc
  samples$thin <- hyper$thin
  samples$n.pars <- hyper$n.pars
  samples$p.names <- hyper$p.names
  out <- suppressMessages(
    summary_recoverone(samples, start, end, ps, digits, prob, verbose)
  )
  return(out)

}

##' Summarise posterior samples
##'
##' This calls severn different variants of summary function to summarise
##' posterior samples
##'
##' @param object posterior samples
##' @param hyper whether to summarise hyper parameters
##' @param start summarise from which iteration.
##' @param end summarise to the end which iteration. For example, set
##' \code{start = 101} and \code{end = 1000}, instructs the function to
##' calculate from 101 to 1000 iteration.
##' @param hmeans a boolean switch indicating to calculate mean of hyper
##' parameters
##' @param hci boolean switch indicating to calculate credible intervals of
##' hyper parameters
##' @param prob a numeric vector, indicating the quantiles to calculate
##' @param recovery a boolean switch indicating if samples are from a recovery
##' study
##' @param ps true parameter values.  This is only for recovery studies
##' @param type calculate type 1 or 2 hyper parameters
##' @param verbose print more information
##' @param digits printing digits
##' @param ... other arguments
##' @export
##' @examples
##' \dontrun{
##' est1 <- summary(hsam[[1]], FALSE)
##' est2 <- summary(hsam[[1]], FALSE, 1, 100)
##' est3 <- summary_one(hsam[[1]], 1, 100, c(.025, .5, .975), verbose = TRUE)
##' est4 <- summary_one(hsam[[1]], 1, 100, c(.025, .5, .975), verbose = F)
##'
##' est5 <- summary_many(hsam, 1, 100, c(.025, .5, .975), FALSE)
##' est6 <- summary_many(hsam, 1, 100, c(.025, .5, .975), TRUE)
##' est7 <- summary(hsam)
##' est8 <- summary(hsam, verbose = TRUE)
##' est9 <- summary(hsam, verbose = FALSE)
##'
##'
##' hest1 <- summary_hyper(hsam, 1, 100, F, F, c(.025, .5, .975), 2, F)
##' hest2 <- summary_hyper(hsam, 1, 100, F, F, c(.05, .5, .9), 2, F)
##' hest3 <- summary(hsam, TRUE)
##' }
summary.model <- function(object, hyper = FALSE, start = 1, end = NA,
  hmeans = FALSE, hci = FALSE, prob = c(0.025, 0.25, 0.5, 0.75, 0.975),
  recovery = FALSE, ps = NA, type = 1, verbose = FALSE, digits = 2, ...) {

  if ( recovery && !is.null(object$theta) ) {
    if (any(is.na(ps))) stop("Some true values are NAs.")
    out <- summary_recoverone(object, start, end, ps, digits, prob, verbose)

  } else if (hyper && recovery)  {
    out <- summary_recoverhyper(object, start, end, ps, type, digits, prob,
      verbose)

  } else if (recovery && is.null(object$theta)) {
    if (any(is.na(ps))) stop("Some true values are NAs.")
    out <- summary_recovermany(object, start, end, ps, digits, prob)

  } else if (hyper) {
    out <- summary_hyper(object, start, end, hmeans, hci, prob, digits, verbose)
  } else if (!is.null(object$theta)) {
    out <- summary_one(object, start, end, prob, verbose)
  } else {
    out <- summary_many(object, start, end, prob, verbose)
  }

  return(out)

}

### Stuck Chains -------------------------------------------------------
##' Which chains get stuck
##'
##' Calculate each chain separately for the mean (across many MCMC iterations)
##' of posterior log-likelihood. If the difference of the means and
##' the median (across chains) of the mean of posterior is greater than the
##' \code{cut}, chains are considered stuck. The default value for \code{cut}
##' is 10. \code{unstick} manually removes stuck chains from posterior samples.
##'
##' @param x posterior samples
##' @param hyper whether x are hierarhcial samples
##' @param cut a criterion deciding if a chain is stuck.
##' @param start start to evaluate from which iteration.
##' @param end end at which iteration for evaeuation.
##' @param verbose a boolean switch to print more information
##' @param digits print how many digits. Default is 2
##' @return \code{PickStuck} gives an index vector; \code{unstick} gives a DMC
##' sample.
##' @examples
##' model <- BuildModel(
##' p.map     = list(A = "1", B = "1", t0 = "1", mean_v = "M", sd_v = "1", st0 = "1"),
##' match.map = list(M = list(s1 = 1, s2 = 2)),
##' factors   = list(S = c("s1", "s2")),
##' constants = c(st0 = 0, sd_v = 1),
##' responses = c("r1", "r2"),
##' type      = "norm")
##' p.vector <- c(A = .75, B = .25, t0 = .2, mean_v.true = 2.5, mean_v.false = 1.5)
##'
##' p.prior <- BuildPrior(
##'   dists = c("tnorm", "tnorm", "beta", "tnorm", "tnorm"),
##'   p1    = c(A = .3, B = .3, t0 = 1, mean_v.true = 1, mean_v.false = 0),
##'   p2    = c(1, 1,   1, 3, 3),
##'   lower = c(0,  0,  0, NA, NA),
##'   upper = c(NA,NA,  1, NA, NA))
##'
##' \dontrun{
##' dat <- simulate(model, 30, ps = p.vector)
##' dmi <- BuildDMI(dat, model)
##' sam <- run(StartNewsamples(5e2, dmi, p.prior))
##' bad <- PickStuck(sam)
##' }
##' @export
PickStuck <- function(x, hyper = FALSE, cut = 10, start = 1,
                      end = NA, verbose = FALSE, digits = 2) {

  if (hyper) {
    out <- PickStuck_hyper(x, cut, start, end, verbose, digits)
  } else if (!is.null(x$theta)) {
    out <- PickStuck_one(x, cut, start, end, verbose, digits)
  } else {
    out <- PickStuck_many(x, cut, start, end)
  }
  return(out)
}


##' @importFrom matrixStats colMeans2
##' @importFrom stats median
PickStuck_hyper <- function(x, cut, start, end, verbose, digits) {

  hyper <- attr(x, "hyper")
  if (is.null(hyper)) stop("Samples are not from a hierarhcial model fit")
  if (is.na(end)) end <- hyper$nmc
  if (end <= start) stop("End must be greater than start")
  iter <- start:end

  pll <- hyper$h_log_likelihoods[iter, ] + hyper$h_summed_log_prior[iter, ]
  mean.ll <- matrixStats::colMeans2(pll)
  names(mean.ll) <- 1:length(mean.ll)
  dev <- -(sort(mean.ll) - median(mean.ll))
  bad <- as.numeric(names(dev)[dev > cut])

  if (verbose) {
    # message("PickStuck_hyper")
    cat("Deviation of mean chain log-likelihood from median of means\n")
    print(round(dev,digits))
    cat("Bad chains: ")
    if (length(bad) == 0) { cat("None\n") } else { cat(bad, "\n") }
  }
  return(bad)
}

##' @importFrom matrixStats colMeans2
##' @importFrom stats median
PickStuck_one <- function(x, cut, start, end, verbose, digits) {

  if (is.na(end)) end <- x$nmc
  if (end <= start) stop("End must be greater than start")
  iter <- start:end

  pll <- x$log_likelihoods[iter, ] + x$summed_log_prior[iter, ]
  mean.ll <- matrixStats::colMeans2(pll)
  names(mean.ll) <- 1:length(mean.ll)
  dev <- -(sort(mean.ll) - median(mean.ll))
  bad <- as.numeric(names(dev)[dev > cut])

  if (verbose) {
    # message("PickStuck_one")
    cat("Deviation of mean chain log-likelihood from median of means\n")
    print(round(dev, digits))
    cat("Bad chains: ")
    if (length(bad) == 0) { cat("None\n") } else { cat(bad, "\n") }
  }
  return(bad)
}

PickStuck_many <- function(x, cut, start, end) {
  sam1 <- x[[1]]
  if (is.na(end)) end <- sam1$nmc
  if (end <= start) stop("End must be greater than start")
  bad <- lapply(x, PickStuck_one, cut, start, end, FALSE)
  return(bad)
}

##' Unstick posterios samples (One subject)
##'
##' @param x posterior samples
##' @param bad a numeric vector, indicating which chains to remove
##' @export
unstick_one <- function(x, bad) {
  cat("unstick_one")

  nchain <- x$n.chains
  if (length(bad) > 0) {
    if (!all(bad %in% 1:nchain)) stop(paste("Index of bad chains must be in 1 to ",
                                            nchain))

    x$theta            <- x$theta[-bad,,]
    x$summed_log_prior <- x$summed_log_prior[,-bad]
    x$log_likelihoods  <- x$log_likelihoods[,-bad]
    x$n.chains         <- x$n.chains - length(bad)
  }

  return(x)
}
