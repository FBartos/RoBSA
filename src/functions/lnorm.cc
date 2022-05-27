#include "lnorm.h"
#include "../source/distributions.h"
#include <util/nainf.h>
#include <JRmath.h>
#include <cmath>

using std::vector;

namespace jags {
  namespace surv {
  
    DLNORM_fun::DLNORM_fun() :ScalarFunction("lnorm_aft_event_lpdf", 3)
    {}
    bool DLNORM_fun::checkParameterValue(vector<double const *> const &args) const
    {
      return(true);
    }
    double DLNORM_fun::evaluate(vector<double const *> const &args) const
    {
			double t      = *args[0];
			double eta    = *args[1];
			double sd     = *args[2];

			return lnorm_aft_log_density(t, eta, sd);
    }

    SLNORM_fun::SLNORM_fun() :ScalarFunction("lnorm_aft_cens_r_lpdf", 3)
    {}
    bool SLNORM_fun::checkParameterValue(vector<double const *> const &args) const
    {
      return(true);
    }
    double SLNORM_fun::evaluate(vector<double const *> const &args) const
    {
			double t      = *args[0];
			double eta    = *args[1];
			double sd     = *args[2];

			return lnorm_aft_log_survival(t, eta, sd);
    }
  }
}
