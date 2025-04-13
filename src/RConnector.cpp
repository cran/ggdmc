#include <ggdmc.hpp>

using namespace Rcpp;

//' Calculate likelihoods
//'
//' These function calculate likelihoods. \code{likelihood_rd} implements
//' the equations in Voss, Rothermund, and Voss (2004). These equations
//' calculate diffusion decision model (Ratcliff & Mckoon, 2008). Specifically,
//' this function implements Voss, Rothermund, and Voss's (2004) equations A1
//' to A4 (page 1217) in C++.
//'
//' @param p_vector_r a parameter vector
//' @param data data model instance
//' @param min_lik minimal likelihood.
//' @return a vector
//' @references Voss, A., Rothermund, K., & Voss, J. (2004).  Interpreting the
//' parameters of the diffusion model: An empirical validation.
//' \emph{Memory & Cognition}, \bold{32(7)}, 1206-1220. \cr\cr
//' Ratcliff, R. (1978). A theory of memory retrival. \emph{Psychological
//' Review}, \bold{85}, 238-255.
//'
//' @examples
//' model <- BuildModel(
//' p.map     = list(A = "1", B = "1", t0 = "1", mean_v = "M", sd_v = "1",
//'             st0 = "1"),
//' match.map = list(M = list(s1 = 1, s2 = 2)),
//' factors   = list(S = c("s1", "s2")),
//' constants = c(st0 = 0, sd_v = 1),
//' responses = c("r1", "r2"),
//' type      = "norm")
//'
//' p.vector <- c(A = .25, B = .35,  t0 = .2, mean_v.true = 1,
//'           mean_v.false = .25)
//' dat <- simulate(model, 1e3,  ps = p.vector)
//' dmi <- BuildDMI(dat, model)
//' den <- likelihood(p.vector, dmi)
//'
//' model <- BuildModel(
//' p.map     = list(a = "1", v = "1", z = "1", d = "1", t0 = "1", sv = "1",
//'             sz = "1", st0 = "1"),
//' constants = c(st0 = 0, d = 0),
//' match.map = list(M = list(s1 = "r1", s2 = "r2")),
//' factors   = list(S = c("s1", "s2")),
//' responses = c("r1", "r2"),
//' type      = "rd")
//'
//' p.vector <- c(a = 1, v = 1, z = 0.5, sz = 0.25, sv = 0.2, t0 = .15)
//' dat <- simulate(model, 1e2, ps = p.vector)
//' dmi <- BuildDMI(dat, model)
//' den <- likelihood (p.vector, dmi)
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector likelihood(Rcpp::NumericVector p_vector_r, List data,
                               double min_lik = 1e-10) {
  // used only by R
  Design *obj0 = new Design(data);
  Likelihood *obj1 = new Likelihood(data, obj0);

  arma::vec p_vector = Rcpp::as<arma::vec>(p_vector_r);
  arma::vec tmp = obj1->likelihood(p_vector);

  Rcpp::NumericVector out(obj0->m_nRT);

  for (size_t i = 0; i < obj0->m_nRT; i++) {
    out[i] = R::fmax2(tmp[i], min_lik);
  }

  delete obj1;
  return out;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix
p_df(Rcpp::NumericVector p_vector_r, std::string cell, std::string mtype,

     std::vector<std::string> pnames, std::vector<std::string> parnames,
     std::vector<std::string> dim0, std::vector<std::string> dim1,
     std::vector<std::string> dim2,

     std::vector<double> allpar, Rcpp::NumericVector model,
     Rcpp::NumericVector isr1, Rcpp::NumericMatrix n1idx, bool n1order) {
  // Used only in random.R

  arma::vec p_vector = Rcpp::as<arma::vec>(p_vector_r);

  arma::ucube arma_model = Rcpp::as<arma::ucube>(model);
  arma::uvec arma_isr1 = Rcpp::as<arma::uvec>(isr1);
  arma::umat arma_n1idx = Rcpp::as<arma::umat>(n1idx);

  Design *obj0 =
      new Design(pnames, parnames, dim0, dim1, dim2, allpar, arma_model);
  Likelihood *obj1 =
      new Likelihood(mtype, arma_isr1, arma_n1idx, n1order, obj0);

  arma::mat pmat = obj1->get_pmat(p_vector, cell);

  Rcpp::NumericMatrix out = Rcpp::wrap(pmat);

  delete obj1; // obj0 is freed in obj1;
  return out;
}

//' Generate Random Responses of the LBA Distribution
//'
//' \code{rlba_norm} used the LBA process to generate response times and
//' responses.
//'
//' @param n is the numbers of observation.
//' @param A start point upper bound, a vector of a scalar.
//' @param b decision threshold, a vector or a scalar.
//' @param mean_v mean drift rate vector
//' @param sd_v standard deviation of drift rate vector
//' @param t0 non-decision time, a vector.
//' @param st0 non-decision time variation, a vector.
//' @param posdrift if exclude negative drift rates
//'
//' @return a n x 2 matrix of response time (first column) and responses (second
//' column).
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix rlba_norm(unsigned int n, const Rcpp::NumericVector A,
                              const Rcpp::NumericVector b,
                              const Rcpp::NumericVector mean_v,
                              const Rcpp::NumericVector sd_v,
                              const Rcpp::NumericVector t0,
                              const Rcpp::NumericVector st0, bool posdrift) {

  // Convert to arma vectors and expand scalars
  auto expand_if_needed = [&mean_v](arma::vec &v) {
    if (v.n_elem == 1)
      v = arma::repmat(v, mean_v.size(), 1);
  };
  arma::vec arma_A = Rcpp::as<arma::vec>(A);
  arma::vec arma_b = Rcpp::as<arma::vec>(b);
  arma::vec arma_t0 = Rcpp::as<arma::vec>(t0);
  arma::vec arma_st0 = Rcpp::as<arma::vec>(st0);
  arma::vec arma_sd_v = Rcpp::as<arma::vec>(sd_v);

  expand_if_needed(arma_A);
  expand_if_needed(arma_b);
  expand_if_needed(arma_t0);
  expand_if_needed(arma_st0);
  expand_if_needed(arma_sd_v);

  // Create LBA object and generate samples
  size_t n_mean_v = mean_v.size();
  std::vector<double> cpp_mean_v(n_mean_v), cpp_sd_v(n_mean_v);

  for (size_t i = 0; i < n_mean_v; ++i) {
    cpp_mean_v[i] = mean_v[i];
    cpp_sd_v[i] = sd_v[i];
  }

  lba obj(arma_A[0], arma_b[0], cpp_mean_v, cpp_sd_v, arma_t0[0], arma_st0[0],
          posdrift);

  arma::mat arma_out(n, 2);
  obj.r(n, arma_out);

  return Rcpp::wrap(arma_out);
}

// [[Rcpp::export]]
NumericMatrix rprior_mat(List prior, unsigned int n) {

  if (n < 1)
    stop("n must be greater or equal to 1");

  Prior *obj = new Prior(prior);
  CharacterVector pnames = prior.attr("names");
  unsigned int npar = pnames.size();

  NumericMatrix out(n, npar);
  for (size_t i = 0; i < n; i++) {
    arma::vec tmp = obj->rprior();

    for (size_t j = 0; j < npar; j++) {
      out(i, j) = tmp[j];
    }
  }

  Rcpp::colnames(out) = pnames;
  return out;
}

// --------------------- Initialization -----------------------------
// [[Rcpp::export]]
List init_mcmc(List data_or_samples, List prior, unsigned int nchain,
               unsigned int nmc, unsigned int thin, unsigned int report,
               double rp, double gammamult, double pm_Hu, double pm_BT,
               bool block, bool add = false, bool is_old = false) {
  // Initialize variables that differ between new/old cases
  List data;
  std::vector<std::string> pnames;
  unsigned int npar;

  // Handle the differences between new and old initialization
  if (is_old) {
    List samples_in = clone(data_or_samples);
    prior = samples_in["p.prior"];
    data = samples_in["data"];
    nchain = samples_in["n.chains"];
    npar = samples_in["n.pars"];
  } else {
    data = data_or_samples;
    npar = prior.size();
  }

  // Common initialization code
  Design *d0 = new Design(data);
  Prior *p0 = new Prior(prior);
  Likelihood *l0 = new Likelihood(data, d0);

  // Initialize Theta differently based on is_old
  Theta *t0 = is_old ? new Theta(data_or_samples, nmc, thin, p0, l0, add)
                     : new Theta(nmc, nchain, npar, thin, p0, l0);

  Sampler *s0 = new Sampler(nchain, npar, gammamult, rp);

  // Get parameter names
  pnames.resize(npar);
  for (size_t i = 0; i < npar; i++) {
    pnames[i] = d0->m_pnames[i];
  }

  // Common sampling loop
  for (size_t i = 1; i < t0->m_nsamp; i++) {
    if (R::runif(0.0, 1.0) < pm_BT) {
      s0->migrate_old(t0);
    } else if (R::runif(0.0, 1.0) < pm_Hu) {
      s0->migrate(t0);
    } else {
      if (block) {
        for (size_t j = 0; j < npar; j++)
          s0->crossover(j, t0);
      } else {
        s0->crossover(t0);
      }
    }
    t0->store(i, report, true);
  }
  Rcout << std::endl;

  // Common output construction
  List out = List::create(
      Named("theta") = t0->m_theta, Named("summed_log_prior") = t0->m_lp,
      Named("log_likelihoods") = t0->m_ll, Named("data") = data,
      Named("p.prior") = prior, Named("start") = t0->m_start_R,
      Named("n.pars") = npar, Named("p.names") = pnames,
      Named("nmc") = is_old ? t0->m_nmc : nmc, Named("thin") = t0->m_thin,
      Named("n.chains") = nchain);

  // Cleanup
  delete t0;
  delete s0;
  return out;
}

// --------------------- Hierarchical versions -----------------------------
static void update_priors(Theta *t, Phi *phi) {
  for (size_t i = 0; i < phi->m_nchain; i++) {
    phi->m_p->m_p0 = phi->m_usephi0.col(i);
    phi->m_p->m_p1 = phi->m_usephi1.col(i);
    t->m_p->m_p0 = phi->m_usephi0.col(i);
    t->m_p->m_p1 = phi->m_usephi1.col(i);
    t->m_uselp[i] = phi->m_p->sumlogprior(t->m_usetheta.col(i));
  }
}

// [[Rcpp::export]]
List init_hier(List data_or_samples, List prior, List lprior, List sprior,
               unsigned int nchain, unsigned int nmc, unsigned int thin,
               unsigned int report, double rp, double gammamult, double pm_Hu,
               double pm_BT, bool block, bool add = false,
               bool is_old = false) {
  // Initialize variables that differ between cases
  List data;
  unsigned int npar;
  unsigned int nsub;

  // Handle the differences between new and old initialization
  if (is_old) {
    List samples_in = clone(data_or_samples);
    List hyper = samples_in.attr("hyper");
    List pprior = hyper["pp.prior"];
    lprior = pprior["location"];
    sprior = pprior["scale"];

    List subject0 = samples_in[0];
    prior = subject0["p.prior"];
    nchain = hyper["n.chains"];
    nsub = samples_in.size();
  } else {
    data = data_or_samples;
    nsub = data.size();
  }

  npar = prior.size();

  // Common initialization code
  Prior *p0 = new Prior(prior);
  Prior *lp = new Prior(lprior);
  Prior *sp = new Prior(sprior);

  std::vector<Design *> ds(nsub);
  std::vector<Prior *> ps(nsub);
  std::vector<Likelihood *> ls(nsub);
  std::vector<Theta *> ts(nsub);

  // Initialize subjects
  for (size_t i = 0; i < nsub; i++) {
    List datai = is_old ? ((List)data_or_samples[i])["data"] : // For old case
                     (List)data_or_samples[i];                 // For new case

    ds[i] = new Design(datai);
    ps[i] = new Prior(prior);
    ls[i] = new Likelihood(datai, ds[i]);

    if (is_old) {
      List subject_data = (List)data_or_samples[i];
      ts[i] = new Theta(subject_data, nmc, thin, ps[i], ls[i], add);
    } else {
      ts[i] = new Theta(nmc, nchain, npar, thin, ps[i], ls[i]);
    }
  }

  // Initialize Phi differently based on is_old
  Phi *phi = is_old ? new Phi(data_or_samples, nmc, nchain, npar, nsub, thin,
                              add, p0, lp, sp)
                    : new Phi(nmc, nchain, npar, nsub, thin, p0, lp, sp, ts);

  Sampler *s0 = new Sampler(nchain, npar, gammamult, rp);
  Rcout << "Start sampling: ";

  // Common sampling loop
  for (size_t i = 1; i < phi->m_nsamp; i++) {
    if (R::runif(0.0, 1.0) < pm_BT) {
      s0->migrate_old(phi, ts);
    } else if (R::runif(0.0, 1.0) < pm_Hu) {
      s0->migrate(phi, ts);
    } else {
      for (size_t j = 0; j < npar; j++)
        s0->crossover(j, phi, ts);
    }

    // Subject-level sampling
    for (size_t k = 0; k < nsub; k++) {
      update_priors(ts[k], phi);

      if (R::runif(0.0, 1.0) < pm_BT) {
        s0->migrate_old(ts[k]);
      } else if (R::runif(0.0, 1.0) < pm_Hu) {
        s0->migrate(ts[k]);
      } else if (block) {
        for (size_t j = 0; j < npar; j++)
          s0->crossover(j, ts[k]);
      } else {
        s0->crossover(ts[k]);
      }

      ts[k]->store(i, report);
    }

    phi->store(i, report);
  }
  Rcout << std::endl;

  // Prepare output
  std::vector<std::string> pnames = prior.attr("names");
  List out(nsub);

  for (size_t i = 0; i < nsub; i++) {
    List datai =
        is_old ? ((List)data_or_samples[i])["data"] : (List)data_or_samples[i];

    out[i] = List::create(
        Named("theta") = ts[i]->m_theta,
        Named("summed_log_prior") = ts[i]->m_lp,
        Named("log_likelihoods") = ts[i]->m_ll, Named("data") = datai,
        Named("p.prior") = prior, Named("start") = ts[i]->m_start_R,
        Named("n.pars") = npar, Named("p.names") = pnames,
        Named("nmc") = is_old ? ts[i]->m_nmc : nmc,
        Named("thin") = ts[i]->m_thin, Named("n.chains") = nchain);
  }

  // Prepare hyperparameters output
  List phi_tmp = List::create(Named("location") = phi->m_phi0,
                              Named("scale") = phi->m_phi1);

  List ppprior =
      List::create(Named("location") = lprior, Named("scale") = sprior);

  List hyper = List::create(
      Named("phi") = phi_tmp, Named("h_summed_log_prior") = phi->m_hlp,
      Named("h_log_likelihoods") = phi->m_hll, Named("pp.prior") = ppprior,
      Named("start") = phi->m_start_R, Named("n.pars") = npar,
      Named("p.names") = pnames, Named("rp") = rp, Named("nmc") = phi->m_nmc,
      Named("thin") = phi->m_thin, Named("n.chains") = nchain);

  out.attr("hyper") = hyper;

  // Cleanup
  delete s0;
  delete phi;
  for (size_t i = 0; i < nsub; i++) {
    delete ts[i];
  }

  return out;
}
