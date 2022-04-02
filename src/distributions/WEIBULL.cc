
#include "WEIBULL.h"
#include <util/nainf.h>
#include <cmath>
#include <rng/RNG.h>
#include <JRmath.h>
#include "distribution_functions.h"

using std::vector;

namespace jags {
	namespace surv {

		DWEIBULL::DWEIBULL() : VectorDist("weibull_aft_event", 2) {}
		SWEIBULL::SWEIBULL() : VectorDist("weibull_aft_cens_r", 2) {}

		bool DWEIBULL::checkParameterLength(vector<unsigned int> const &len) const
		{
			return true;
		}
		bool SWEIBULL::checkParameterLength(vector<unsigned int> const &len) const
		{
			return true;
		}

		bool DWEIBULL::checkParameterValue(vector<double const *> const &par,
						vector<unsigned int> const &len) const
		{
			return *par[1] > 0.0;
		}
		bool SWEIBULL::checkParameterValue(vector<double const *> const &par,
						vector<unsigned int> const &len) const
		{
			return *par[1] > 0.0;
		}

		double DWEIBULL::logDensity(double const *x, unsigned int length, PDFType type,
					vector<double const *> const &par,
					vector<unsigned int> const &len,
					double const *lower, double const *upper) const
		{
			double t      = *x;
			double eta    = *par[0];
			double shape  = *par[1];

			return weibull_aft_log_density(t, eta, shape);
		}
		double SWEIBULL::logDensity(double const *x, unsigned int length, PDFType type,
					vector<double const *> const &par,
					vector<unsigned int> const &len,
					double const *lower, double const *upper) const
		{
			double t      = *x;
			double eta    = *par[0];
			double shape  = *par[1];

			return weibull_aft_log_survival(t, eta, shape);
		}

		void DWEIBULL::randomSample(double *x, unsigned int length,
					vector<double const *> const &par,
					vector<unsigned int> const &len,
					double const *lower, double const *upper,
					RNG *rng) const
		{
			// not implemented
		}
		void SWEIBULL::randomSample(double *x, unsigned int length,
					vector<double const *> const &par,
					vector<unsigned int> const &len,
					double const *lower, double const *upper,
					RNG *rng) const
		{
			// not implemented
		}

		void DWEIBULL::support(double *lower, double *upper, unsigned int length,
				vector<double const *> const &par,
				vector<unsigned int> const &len) const
		{
			for (unsigned int i = 0; i < length; ++i) {
				lower[i] = 0;
				upper[i] = JAGS_POSINF;
			}
		}
		void SWEIBULL::support(double *lower, double *upper, unsigned int length,
				vector<double const *> const &par,
				vector<unsigned int> const &len) const
		{
			for (unsigned int i = 0; i < length; ++i) {
				lower[i] = 0;
				upper[i] = JAGS_POSINF;
			}
		}

		unsigned int DWEIBULL::length(vector<unsigned int> const &len) const
		{
			return 1;
		}
		unsigned int SWEIBULL::length(vector<unsigned int> const &len) const
		{
			return 1;
		}


		void DWEIBULL::typicalValue(double *x, unsigned int length,
					vector<double const *> const &par,
					vector<unsigned int> const &len,
					double const *lower, double const *upper) const
		{
			// not implemented
		}
		void SWEIBULL::typicalValue(double *x, unsigned int length,
					vector<double const *> const &par,
					vector<unsigned int> const &len,
					double const *lower, double const *upper) const
		{
			// not implemented
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

