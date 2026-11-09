#include "WEIBULL.h"
#include <util/nainf.h>
#include <cmath>
#include <rng/RNG.h>
#include <JRmath.h>
#include "../source/distributions.h"

using std::vector;

namespace jags {
	namespace surv {

		DWEIBULL::DWEIBULL() : VectorDist("weibull_aft_event", 2) {}
		SWEIBULL::SWEIBULL() : VectorDist("weibull_aft_cens_r", 2) {}

		bool DWEIBULL::checkParameterLength(vector<unsigned long> const &len) const
		{
			return true;
		}
		bool SWEIBULL::checkParameterLength(vector<unsigned long> const &len) const
		{
			return true;
		}

		bool DWEIBULL::checkParameterValue(vector<double const *> const &par,
						vector<unsigned long> const &len) const
		{
			return *par[1] > 0.0;
		}
		bool SWEIBULL::checkParameterValue(vector<double const *> const &par,
						vector<unsigned long> const &len) const
		{
			return *par[1] > 0.0;
		}

		double DWEIBULL::logDensity(double const *x, PDFType type,
					vector<double const *> const &par,
					vector<unsigned long> const &len) const
		{
			double t      = *x;
			double eta    = *par[0];
			double shape  = *par[1];

			return weibull_aft_log_density(t, eta, shape);
		}
		double SWEIBULL::logDensity(double const *x,  PDFType type,
					vector<double const *> const &par,
					vector<unsigned long> const &len) const
		{
			double t      = *x;
			double eta    = *par[0];
			double shape  = *par[1];

			return weibull_aft_log_survival(t, eta, shape);
		}

		void DWEIBULL::randomSample(double *x, 
					vector<double const *> const &par,
					vector<unsigned long> const &len,
					RNG *rng) const
		{
			// not implemented
		}
	    
		void SWEIBULL::randomSample(double *x, 
					vector<double const *> const &par,
					vector<unsigned long> const &len,
					RNG *rng) const
		{
			// not implemented
		}

		void DWEIBULL::support(double *lower, double *upper, 
				vector<double const *> const &par,
				vector<unsigned long> const &len) const
		{
		    *lower = 0;
		    *upper = JAGS_POSINF;
		}
	    
		void SWEIBULL::support(double *lower, double *upper, 
				vector<double const *> const &par,
				vector<unsigned long> const &len) const
		{
		    *lower = 0;
		    *upper = JAGS_POSINF;
		}

		unsigned long DWEIBULL::length(vector<unsigned long> const &len) const
		{
			return 1;
		}
		unsigned long SWEIBULL::length(vector<unsigned long> const &len) const
		{
			return 1;
		}

		bool DWEIBULL::isSupportFixed(vector<bool> const &fixmask) const
		{
			return true;
		}
		bool SWEIBULL::isSupportFixed(vector<bool> const &fixmask) const
		{
			return true;
		}
	}
}

