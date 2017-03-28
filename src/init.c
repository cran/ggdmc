#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

extern SEXP ggdmc_getAccumulatorMatrix(SEXP, SEXP, SEXP, SEXP);
extern SEXP ggdmc_ddmc(SEXP, SEXP, SEXP, SEXP);
extern SEXP ggdmc_ddmc_parallel(SEXP, SEXP, SEXP, SEXP);
extern SEXP ggdmc_g_minus(SEXP);
extern SEXP ggdmc_g_plus(SEXP);
extern SEXP ggdmc_g_minus_parallel(SEXP);
extern SEXP ggdmc_g_plus_parallel(SEXP);
extern SEXP ggdmc_initialise_data(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP ggdmc_initialise_hyper(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP ggdmc_dprior(SEXP, SEXP);
extern SEXP ggdmc_rprior(SEXP, SEXP);
extern SEXP ggdmc_run_data(SEXP, SEXP, SEXP);
extern SEXP ggdmc_run_data_parallel(SEXP, SEXP, SEXP);
extern SEXP ggdmc_run_hyper(SEXP, SEXP);
extern SEXP ggdmc_run_hyper_parallel(SEXP, SEXP);
extern SEXP ggdmc_assign_pp_pLists(SEXP, SEXP);
extern SEXP ggdmc_summed_log_likelihood(SEXP, SEXP);
extern SEXP ggdmc_summed_log_likelihood_parallel(SEXP, SEXP);
extern SEXP ggdmc_summed_log_prior(SEXP, SEXP);
extern SEXP rtn_wrapper(SEXP, SEXP, SEXP, SEXP);
extern SEXP dtn_wrapper(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"ggdmc_getAccumulatorMatrix", (DL_FUNC) &ggdmc_getAccumulatorMatrix, 4},
  {"ggdmc_ddmc", (DL_FUNC) &ggdmc_ddmc, 4},
  {"ggdmc_ddmc_parallel", (DL_FUNC) &ggdmc_ddmc_parallel, 4},
  {"ggdmc_g_minus", (DL_FUNC) &ggdmc_g_minus, 1},
  {"ggdmc_g_plus", (DL_FUNC) &ggdmc_g_minus, 1},
  {"ggdmc_g_minus_parallel", (DL_FUNC) &ggdmc_g_minus_parallel, 1},
  {"ggdmc_g_plus_parallel", (DL_FUNC) &ggdmc_g_plus_parallel, 1},
  {"ggdmc_initialise_data", (DL_FUNC) &ggdmc_initialise_data, 7},
  {"ggdmc_initialise_hyper", (DL_FUNC) &ggdmc_initialise_hyper, 10},
  {"ggdmc_dprior", (DL_FUNC) &ggdmc_dprior, 2},
  {"ggdmc_rprior", (DL_FUNC) &ggdmc_rprior, 2},
  {"ggdmc_run_data", (DL_FUNC) &ggdmc_run_data, 3},
  {"ggdmc_run_data_parallel", (DL_FUNC) &ggdmc_run_data_parallel, 3},
  {"ggdmc_run_hyper", (DL_FUNC) &ggdmc_run_hyper, 2},
  {"ggdmc_run_hyper_parallel", (DL_FUNC) &ggdmc_run_hyper, 2},
  {"ggdmc_assign_pp_pLists", (DL_FUNC) &ggdmc_assign_pp_pLists, 2},
  {"ggdmc_summed_log_likelihood", (DL_FUNC) &ggdmc_summed_log_likelihood, 2},
  {"ggdmc_summed_log_likelihood_parallel", (DL_FUNC) &ggdmc_summed_log_likelihood_parallel, 2},
  {"ggdmc_summed_log_prior", (DL_FUNC) &ggdmc_summed_log_prior, 2},
  {"rtn_wrapper", (DL_FUNC) &rtn_wrapper, 4},
  {"dtn_wrapper", (DL_FUNC) &dtn_wrapper, 6},
  {NULL, NULL, 0}
};

void R_init_ggdmc(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}

