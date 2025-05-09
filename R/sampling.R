### Sampling -------------------------------------------------------------------
run_one <- function (data, prior, nchain, nmc, thin, report, rp, gammamult,
                     pm_Hu, pm_BT, block)
{
  out <- init_mcmc(data, prior, nchain, nmc, thin, report, rp, gammamult, pm_Hu,
                  pm_BT, block)
  pnames <- GetPNames(attr(out$data, "model"))
  dimnames(out$theta) <- list(pnames, NULL, NULL)
  return(out)
}


run_hier <- function (prior, lprior, sprior, data, nchain, nmc, thin, report,
                      rp, gammamult, pm_Hu, pm_BT, block, add, is_old)
{
  out <- init_hier(data, prior, lprior, sprior, nchain, nmc, thin, report,
                   rp, gammamult, pm_Hu, pm_BT, block, add, is_old)

  pnames <- GetPNames(attr(out[[1]]$data, "model"))
  phi1_tmp <- attr(out, "hyper")$phi[[1]]
  phi2_tmp <- attr(out, "hyper")$phi[[2]]
  dimnames(phi1_tmp) <- list(pnames, NULL, NULL)
  dimnames(phi2_tmp) <- list(pnames, NULL, NULL)
  attr(out, "hyper")$phi[[1]] <- phi1_tmp
  attr(out, "hyper")$phi[[2]] <- phi2_tmp

  for(i in 1:length(out))
  {
    dimnames(out[[i]]$theta) <- list(pnames, NULL, NULL)
    class(out[[i]]) <- c("model", "list")
  }

  names(out) <- names(data)

  return(out)
}


##' @importFrom parallel parLapply
##' @importFrom parallel makeCluster
##' @importFrom parallel mclapply
##' @importFrom parallel stopCluster
run_many <- function(data, prior, nchain, nmc, thin, report, rp, gammamult,
                     pm0, pm_BT, block, ncore)
{

  if (get_os() == "windows" & ncore > 1)
  {
    cl  <- parallel::makeCluster(ncore)
    message("Run subjects in parallel")
    out <- parallel::parLapply(cl = cl, X = data,
                               init_mcmc, prior, nchain, nmc, thin, report,
                               rp, gammamult, pm0, pm_BT, block)
    parallel::stopCluster(cl)
  }
  else if (ncore > 1)
  {
    message("Run subjects in parallel")
    out <- parallel::mclapply(data, init_mcmc, prior, nchain, nmc, thin, report,
                              rp, gammamult, pm0, pm_BT, block, mc.cores = ncore)
  }
  else
  {
    message("Run subjects with lapply")
    out <- lapply(data, init_mcmc, prior, nchain, nmc, thin, report, rp,
                  gammamult, pm0, pm_BT, block)
  }

  for(i in 1:length(out))
    dimnames(out[[i]]$theta) <- list(out[[i]]$p.names, NULL, NULL)

  return(out)
}

##' @importFrom parallel parLapply
##' @importFrom parallel makeCluster
##' @importFrom parallel mclapply
##' @importFrom parallel stopCluster
rerun_many <- function(samples, nmc, thin, report, rp, gammamult, pm_Hu, pm_BT,
                       block, add, ncore, is_old)
{

  if (get_os() == "windows" & ncore > 1)
  {
    cl  <- parallel::makeCluster(ncore)
    message("Run subjects in parallel (parLapply)")
    n_subject <- length(samples)
    out <- parallel::parLapply(cl = cl, X = seq_len(n_subject), function(i) {
      samples_in <- samples[[i]]
      p_prior <- samples[[i]]$p.prior
      n_chain <- samples[[i]]$n.chains
      init_mcmc(samples[[i]], p_prior, n_chain, nmc, thin, report, rp,
                gammamult, pm_Hu, pm_BT, block, add, is_old)
    })

    
    parallel::stopCluster(cl)
  }
  else if (ncore > 1)
  {
    message("Run subjects in parallel (mclapply)")
    n_subject <- length(samples)
    out <- parallel::mclapply(seq_len(n_subject), function(i) {
      samples_in <- samples[[i]]
      p_prior <- samples[[i]]$p.prior
      n_chain <- samples[[i]]$n.chains
      init_mcmc(samples[[i]], p_prior, n_chain, nmc, thin, report, rp,
                gammamult, pm_Hu, pm_BT, block, add, is_old)
    }, mc.cores = ncore)

  }
  else
  {
    message("Run subjects with lapply")
    n_subject <- length(samples)
    out <- lapply(seq_len(n_subject), function(i) {
      samples_in <- samples[[i]]
      p_prior <- samples[[i]]$p.prior
      n_chain <- samples[[i]]$n.chains
      init_mcmc(samples[[i]], p_prior, n_chain, nmc, thin, report, rp,
                gammamult, pm_Hu, pm_BT, block, add, is_old)
      })
  }

  for(i in 1:length(out))
    dimnames(out[[i]]$theta) <- list(out[[i]]$p.names, NULL, NULL)

  return(out)
}

##' Start new model fits
##'
##' Fit a hierarchical or a fixed-effect model, using Bayeisan
##' optimisation.  We use a specific type of pMCMC algorithm, the DE-MCMC. This
##' particular sampling method includes crossover and two different migration
##' operators. The migration operators are similar to random-walk algorithm.
##' They wouold be less efficient to find the target parameter space, if been
##' used alone.
##'
##' @param data data model instance(s)
##' @param samples posterior samples.
##' @param prior prior objects.  For hierarchical model, this must be a
##' list with three sets of prior distributions. Each is respectively named,
##' "pprior", "location", and "scale".
##' @param nmc number of Monte Carlo samples
##' @param thin thinning length
##' @param nchain number of chains
##' @param report progress report interval
##' @param rp tuning parameter 1
##' @param gammamult tuning parameter 2. This is the step size.
##' @param pm_Hu probability of migration type 0 (Hu & Tsui, 2010)
##' @param pm_BT probability of migration type 1 (Turner et al., 2013)
##' @param block Only for hierarchical modeling. A Boolean switch for update one
##' parameter at a time
##' @param ncore Only for non-hierarchical, fixed-effect models with many
##' subjects.
##' @param add Boolean whether to add new samples
##' @param is_old start sampling from a DMI or fit samples
##'
##' @export
StartNewsamples <- function(data, prior=NULL, nmc=2e2, thin=1, nchain=NULL,
                            report=1e2, rp=.001, gammamult=2.38, pm_Hu=.05,
                            pm_BT=.05, block=TRUE, ncore=1, add = FALSE, is_old = FALSE)
{
  if ( !is.null(prior) &&
       all(c("pprior", "location", "scale") %in% names(prior)) &&
       length(prior) == 3 )
  {
    nchain <- CheckHyperDMI(data, prior, nchain)
    checklba(data[[1]])


    message("Hierarchical model")
    out <- run_hier(prior[[1]], prior[[2]], prior[[3]], data, nchain, nmc, thin,
                    report, rp, gammamult, pm_Hu, pm_BT, block, add, is_old)
  }
  else if ( is.data.frame(data) )
  {
    nchain <- CheckDMI(data, prior, nchain)
    checklba(data)

    message("Non-hierarchical model")
    out <- run_one(data, prior, nchain, nmc, thin, report, rp, gammamult, pm_Hu,
                   pm_BT, block)
  }
  else
  {
    for (i in 1:length(data)) nchain <- CheckDMI(data[[i]], prior, nchain)
    checklba(data[[1]])

    message("Non-hierarchical model fitting many subjects")
    out <- run_many(data, prior, nchain, nmc, thin, report, rp, gammamult, pm_Hu,
                    pm_BT, block, ncore)
  }

  class(out) <- c("model", "list")
  cat("\n")
  return(out)
}

##' @rdname StartNewsamples
##' @export
run <-  function(samples, nmc=5e2, thin=1, report=1e2, rp=.001,
                 gammamult=2.38, pm_Hu=0, pm_BT=0, block=TRUE, ncore=1,
                 add=FALSE, is_old = TRUE)
{
  hyper <- attr(samples, "hyper")

  if ( !is.null(hyper) )
  {
    subject1 <- samples[[1]]
    pprior <- subject1$p.prior
    lprior <- hyper$pp.prior$location
    sprior <- hyper$pp.prior$scale
    nchain <- subject1$n.chains
    
    out <- init_hier(samples, pprior, lprior, sprior, nchain, nmc, thin, report, rp, gammamult, pm_Hu, pm_BT,
                        block, add, is_old)

    pnames   <- GetPNames(attr(out[[1]]$data, "model"))
    phi1_tmp <- attr(out, "hyper")$phi[[1]]
    phi2_tmp <- attr(out, "hyper")$phi[[2]]
    dimnames(phi1_tmp) <- list(pnames, NULL, NULL)
    dimnames(phi2_tmp) <- list(pnames, NULL, NULL)
    attr(out, "hyper")$phi[[1]] <- phi1_tmp
    attr(out, "hyper")$phi[[2]] <- phi2_tmp

    for(i in 1:length(out))
    {
      dimnames(out[[i]]$theta) <- list(pnames, NULL, NULL)
      class(out[[i]]) <- c("model", "list")
    }

    names(out) <- names(samples)
  }
  else if (any(names(samples) == "theta"))
  {
    out <- init_mcmc(samples, samples$p.prior, samples$n.chains, nmc, thin, report, rp, gammamult, pm_Hu, pm_BT,
                    block, add, is_old)
    pnames <- GetPNames(attr(out$data, "model"))
    dimnames(out$theta) <- list(pnames, NULL, NULL)
  }
  else
  {
    out <- rerun_many(samples, nmc, thin, report, rp, gammamult, pm_Hu, pm_BT,
                      block, add, ncore, is_old)
  }

  class(out) <- c("model", "list")
  cat("\n")
  return(out)
}

## Utilities ------------------------------------------------------------------
CheckDMI <- function(data = NULL, prior = NULL, nchain = NULL)
{
  if (is.null(data))        stop("No data model instance")
  if (!is.data.frame(data)) stop("Data must be a data frame")

  model <- attr(data, "model")
  npar  <- length(GetPNames(model))
  if (is.null(nchain)) nchain <- 3*npar

  if (is.null(model)) stop("Must specify a model")
  if (is.null(prior)) stop("Must specify prior distributions")
  if (length(prior) != npar)
  {
    message("Must names pprior, location and scale priors")
    stop("data prior incompatible with model")
  }
  return(nchain)
}

CheckHyperDMI <- function(data = NULL, prior = NULL, nchain = NULL)
{
  ## data
  if (is.null(data)) stop("No data")
  if (!is.list(data)) stop ("data must be a list")
  if (is.data.frame(data)) stop("data must not be a data.frame")

  ## Prior
  for (i in 1:3) if (is.null(prior[[i]])) stop("No prior distribution")
  if ( length(prior[[1]]) < length(prior[[2]]) )
    stop("Location and scale priors differ in the numbers of parameters")

  model1 <- attr(data[[1]], "model")
  pnames <- GetPNames(model1)
  isppriorok <- pnames %in% names(prior[[1]])
  islpriorok <- pnames %in% names(prior[[2]])
  isspriorok <- pnames %in% names(prior[[3]])
  if (!all(isppriorok)) stop("data prior incompatible with model")
  if (!all(islpriorok)) stop("location prior incompatible with model")
  if (!all(isppriorok)) stop("scale prior incompatible with model")

  ## nchain
  if (is.null(nchain)) {
    nchain <- 3*length(pnames)
    message("Use default ", nchain, " chains")
  }

  return(nchain)
}
