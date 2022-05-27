#include "weibull.h"
#include "../source/distributions.h"
#include <util/nainf.h>
#include <JRmath.h>
#include <cmath>

using std::vector;

namespace jags {
  namespace surv {
  
    DWEIBULL_fun::DWEIBULL_fun() :ScalarFunction("weibull_aft_event_lpdf", 3)
    {}
    bool DWEIBULL_fun::checkParameterValue(vector<double const *> const &args) const
    {
      return(true);
    }
    double DWEIBULL_fun::evaluate(vector<double const *> const &args) const
    {
			double t      = *args[0];
			double eta    = *args[1];
			double shape  = *args[2];

			return weibull_aft_log_density(t, eta, shape);
    }

    SWEIBULL_fun::SWEIBULL_fun() :ScalarFunction("weibull_aft_cens_r_lpdf", 3)
    {}
    bool SWEIBULL_fun::checkParameterValue(vector<double const *> const &args) const
    {
      return(true);
    }
    double SWEIBULL_fun::evaluate(vector<double const *> const &args) const
    {
			double t      = *args[0];
			double eta    = *args[1];
			double shape  = *args[2];

			return weibull_aft_log_survival(t, eta, shape);
    }
  }
}
