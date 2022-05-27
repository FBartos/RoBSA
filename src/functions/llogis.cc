#include "llogis.h"
#include "../source/distributions.h"
#include <util/nainf.h>
#include <JRmath.h>
#include <cmath>

using std::vector;

namespace jags {
  namespace surv {
  
    DLLOGIS_fun::DLLOGIS_fun() :ScalarFunction("llogis_aft_event_lpdf", 3)
    {}
    bool DLLOGIS_fun::checkParameterValue(vector<double const *> const &args) const
    {
      return(true);
    }
    double DLLOGIS_fun::evaluate(vector<double const *> const &args) const
    {
			double t      = *args[0];
			double eta    = *args[1];
			double shape  = *args[2];

			return llogis_aft_log_density(t, eta, shape);
    }

    SLLOGIS_fun::SLLOGIS_fun() :ScalarFunction("llogis_aft_cens_r_lpdf", 3)
    {}
    bool SLLOGIS_fun::checkParameterValue(vector<double const *> const &args) const
    {
      return(true);
    }
    double SLLOGIS_fun::evaluate(vector<double const *> const &args) const
    {
			double t      = *args[0];
			double eta    = *args[1];
			double shape  = *args[2];

			return llogis_aft_log_survival(t, eta, shape);
    }
  }
}
