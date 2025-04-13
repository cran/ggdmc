#include <ggdmc.hpp>

using namespace Rcpp;

Prior::Prior (List & pprior)
{
  std::vector<std::string> pnames = pprior.attr("names");
  m_npar = pnames.size();

  arma::vec p0(m_npar), p1(m_npar), l(m_npar), u(m_npar);
  arma::uvec d(m_npar), lg(m_npar);

  for (size_t i = 0; i < m_npar; i++) {
    List a_list = pprior[pnames[i]];
    unsigned int a_dist = a_list.attr("dist");

    d[i]  = a_dist;
    p0[i] = a_list[0];
    p1[i] = a_list[1];
    l[i]  = a_list[2];
    u[i]  = a_list[3];
    lg[i] = a_list[4];
  }
  m_d  = d;
  m_p0 = p0;
  m_p1 = p1;
  m_l  = l;
  m_u  = u;
  m_lg = lg;
}

Prior::~Prior() {}

void Prior::dprior(double * pvector, double * out)
{
  double x, l, u;

  for (size_t i = 0; i < m_npar; i++)
  {
    // NA go here; NA will be converted to 0 (unsigned int type)
    if ( ISNAN(m_p1[i]) || ISNAN(m_d[i]) ) {
      out[i] = m_lg[i] ? R_NegInf : 0;
    } else if ( m_d[i] == TNORM ) {

      l = ISNAN(m_l[i]) ? R_NegInf : m_l[i];
      u = ISNAN(m_u[i]) ? R_PosInf : m_u[i];

      tnorm * obj = new tnorm(m_p0[i], m_p1[i], l, u, m_lg[i]);
      out[i] = obj->d(pvector[i]);
      delete obj;

    } else if ( m_d[i] == BETA_LU ) {

      l = ISNAN(m_l[i]) ? 0 : m_l[i]; // In case the user enters NAs.
      u = ISNAN(m_u[i]) ? 1 : m_u[i];

      x = (pvector[i] - l) / (u -  l);

      // Note m_l differs from m_lg !!!
      out[i] = !m_lg[i] ? R::dbeta(x, m_p0[i], m_p1[i], false) / (u - l) :
                          R::dbeta(x, m_p0[i], m_p1[i], true) - std::log(u - l);

    } else if ( m_d[i] == GAMMA_L ) {

      l = ISNAN(m_l[i]) ? 0 : m_l[i];
      x = ( !R_FINITE(l) ) ? pvector[i] : pvector[i] - l;
      out[i] = R::dgamma(x, m_p0[i], m_p1[i], m_lg[i]);

    } else if ( m_d[i] == LNORM_L ) {

      l = ISNAN(m_l[i]) ? 0 : m_l[i];
      x = ( !R_FINITE(l) ) ? pvector[i] : pvector[i] - l;
      out[i] = R::dlnorm(x, m_p0[i], m_p1[i], m_lg[i]);

    } else if ( m_d[i] == 5 ) {

      out[i] = R::dunif(pvector[i], m_p0[i], m_p1[i], m_lg[i]);

    } else if ( m_d[i] == CONSTANT ) {

      out[i] = m_lg[i] ? R_NegInf : 0;

    } else if ( m_d[i] == TNORM_TAU ) {
      l = ISNAN(m_l[i]) ? R_NegInf : m_l[i];
      u = ISNAN(m_u[i]) ? R_PosInf : m_u[i];

      tnorm * obj = new tnorm(m_p0[i], m_p1[i], l, u, m_lg[i]);
      out[i] = obj->d2(pvector[i]);
      delete obj;

    } else {
      Rcpp::Rcout << "Distribution type undefined" << "\n";
      out[i] = m_lg[i] ? R_NegInf : 0;
    }
  }

}

arma::vec Prior::dprior(arma::vec pvector)
{
  double * pvec = new double[m_npar];
  double * tmp  = new double[m_npar];

  for (size_t i = 0; i < m_npar; i++) pvec[i] = pvector[i];

  dprior(pvec, tmp);

  arma::vec out(m_npar);
  for (size_t i = 0; i < m_npar; i++)
  {
    if ( !R_FINITE(tmp[i]) )
    {
      out[i] = m_lg[i] ? -23.02585 : 1e-10; // critical to hierarchical
    }
    else
    {
      out[i] = tmp[i];
    }
  }

  delete [] pvec;
  delete [] tmp;

  return(out);

}

arma::vec Prior::rprior()
// Used by ininitlise & R's rprior
{
  // replace DMC modified r-function; used in initialise.cpp internally
  double l, u;
  arma::vec out(m_npar); out.fill(NA_REAL);

  // [p1 p2]: [mean sd]; [shape1 shape2]; [shape scale]; [meanlog sdlog]
  for (size_t i = 0; i < m_npar;  i++) {
    if ( ISNAN(m_d[i]) ) {
      out[i] = NA_REAL;

    } else if ( m_d[i] == 1 ) {         // tnorm
      l = ISNAN(m_l[i]) ? R_NegInf : m_l[i];
      u = ISNAN(m_u[i]) ? R_PosInf : m_u[i];

      tnorm * obj = new tnorm(m_p0[i], m_p1[i], l, u);
      out[i] = obj->r();
      delete obj;

    } else if ( m_d[i] == 2 ) {  // beta_ul
      l = ISNAN(m_l[i]) ? 0 : m_l[i];
      u = ISNAN(m_u[i]) ? 1 : m_u[i];
      out[i] = l + R::rbeta(m_p0[i], m_p1[i]) * (u - l);


    } else if ( m_d[i] == 3 ) {  // gamma_l
      l = ISNAN(m_l[i]) ? 0 : m_l[i];
      out[i] = R::rgamma(m_p0[i], m_p1[i]) + l;

    } else if ( m_d[i] == 4 ) {  // lnorm_l
      l = ISNAN(m_l[i]) ? 0 : m_l[i];
      out[i] = R::rlnorm(m_p0[i], m_p1[i]) + l;

    } else if ( m_d[i] == 5 ) {
      out[i] = R::runif(m_p0[i], m_p1[i]);
    } else if ( m_d[i] == 6 ){  // constant
      out[i] = m_p0[i];

    } else if ( m_d[i] == 7 ) { // tnorm2

      Rcout << "Distribution type not supported\n";

    } else {
      Rcout << "Distribution type not supported\n";
      out[i] = NA_REAL;
    }
  }

  return out;
}

double Prior::sumlogprior(arma::vec pvector)
{
  arma::vec out = dprior(pvector);
  return arma::accu(out);
}
