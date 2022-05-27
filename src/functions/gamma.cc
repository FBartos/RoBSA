#include "gamma.h"
#include "../source/distributions.h"
#include <util/nainf.h>
#include <JRmath.h>
#include <cmath>

using std::vector;

namespace jags {
  namespace surv {
  
    DGAMMA_fun::DGAMMA_fun() :ScalarFunction("gamma_aft_event_lpdf", 3)
    {}
    bool DGAMMA_fun::checkParameterValue(vector<double const *> const &args) const
    {
      return(true);
    }
    double DGAMMA_fun::evaluate(vector<double const *> const &args) const
    {
			double t      = *args[0];
			double eta    = *args[1];
			double shape  = *args[2];

			return gamma_aft_log_density(t, eta, shape);
    }

    SGAMMA_fun::SGAMMA_fun() :ScalarFunction("gamma_aft_cens_r_lpdf", 3)
    {}
    bool SGAMMA_fun::checkParameterValue(vector<double const *> const &args) const
    {
      return(true);
    }
    double SGAMMA_fun::evaluate(vector<double const *> const &args) const
    {
			double t      = *args[0];
			double eta    = *args[1];
			double shape  = *args[2];

			return gamma_aft_log_survival(t, eta, shape);
    }
  }
}
