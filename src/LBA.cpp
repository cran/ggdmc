#include <ggdmc.hpp>

using namespace Rcpp;

arma::vec fptpdf(arma::vec rt, double A, double b, double mean_v, double sd_v,
                 double t0, double st0, bool posdrift)
{
  lba * obj = new lba(A, b, mean_v, sd_v, t0, st0, posdrift, rt);
  arma::vec out(obj->m_nrt);

  if(!obj->ValidateParams(false))
  {
    out.fill(1e-10);
  }
  else
  {
    out = obj->d();
  }

  delete obj;
  return out;
}

arma::vec fptcdf(arma::vec rt, double A, double b, double mean_v, double sd_v,
                 double t0, double st0, bool posdrift)
{
  lba * obj = new lba(A, b, mean_v, sd_v, t0, st0, posdrift, rt);
  arma::vec out(obj->m_nrt);

  if(!obj->ValidateParams(false))
  {
    out.fill(1e-10);
  }
  else
  {
    out = obj->p();
  }

  delete obj;
  return out;
}

arma::vec n1PDFfixedt0(arma::vec rt, arma::vec A, arma::vec b, arma::vec mean_v,
                       arma::vec sd_v, arma::vec t0, arma::vec st0,
                       bool posdrift) {

  unsigned int nmean_v = mean_v.n_elem;  // Number of accumulators/responses.
  unsigned int n       = rt.n_elem;      // Number of trials
  unsigned int nsd_v   = sd_v.n_elem;    // Check for matrix operations
  unsigned int nA      = A.n_elem;
  unsigned int nb      = b.n_elem;
  unsigned int nt0     = t0.n_elem;
  unsigned int nst0    = st0.n_elem;      // reduntant

  if (nsd_v == 1) sd_v = arma::repmat(sd_v, nmean_v, 1);
  if (nA    == 1) A    = arma::repmat(A,    nmean_v, 1);
  if (nb    == 1) b    = arma::repmat(b,    nmean_v, 1);
  if (nt0   == 1) t0   = arma::repmat(t0,   nmean_v, 1);
  if (nst0  == 1) st0  = arma::repmat(st0,  nmean_v, 1);

  arma::vec onevec = arma::ones<arma::vec>(n);
  arma::vec node1den = fptpdf(rt, A[0], b[0], mean_v[0], sd_v[0], t0[0], st0[0],
                              posdrift);

  if (nmean_v > 1)
  {
    for (size_t i = 1; i < nmean_v; i++)
    {
      node1den = node1den % (onevec - fptcdf(rt, A[i], b[i], mean_v[i],
                                             sd_v[i], t0[i], st0[i], posdrift));
    }
  }

  return node1den;
}

