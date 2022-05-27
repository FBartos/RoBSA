#include "exp.h"
#include "../source/distributions.h"
#include <util/nainf.h>
#include <JRmath.h>
#include <cmath>

using std::vector;

namespace jags {
  namespace surv {
  
    DEXP_fun::DEXP_fun() :ScalarFunction("exp_aft_event_lpdf", 2)
    {}
    bool DEXP_fun::checkParameterValue(vector<double const *> const &args) const
    {
      return(true);
    }
    double DEXP_fun::evaluate(vector<double const *> const &args) const
    {
			double t      = *args[0];
			double eta    = *args[1];

			return exp_aft_log_density(t, eta);
    }

    SEXP_fun::SEXP_fun() :ScalarFunction("exp_aft_cens_r_lpdf", 2)
    {}
    bool SEXP_fun::checkParameterValue(vector<double const *> const &args) const
    {
      return(true);
    }
    double SEXP_fun::evaluate(vector<double const *> const &args) const
    {
			double t      = *args[0];
			double eta    = *args[1];

			return exp_aft_log_survival(t, eta);
    }
  }
}
