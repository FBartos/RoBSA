
#include "GAMMA.h"
#include <util/nainf.h>
#include <cmath>
#include <rng/RNG.h>
#include <JRmath.h>
#include "distribution_functions.h"

using std::vector;

namespace jags {
	namespace surv {

		DGAMMA::DGAMMA() : VectorDist("gamma_aft_event", 2) {}
		SGAMMA::SGAMMA() : VectorDist("gamma_aft_cens_r", 2) {}

		bool DGAMMA::checkParameterLength(vector<unsigned int> const &len) const
		{
			return true;
		}
		bool SGAMMA::checkParameterLength(vector<unsigned int> const &len) const
		{
			return true;
		}

		bool DGAMMA::checkParameterValue(vector<double const *> const &par,
						vector<unsigned int> const &len) const
		{
			return *par[1] > 0.0;
		}
		bool SGAMMA::checkParameterValue(vector<double const *> const &par,
						vector<unsigned int> const &len) const
		{
			return *par[1] > 0.0;
		}

		double DGAMMA::logDensity(double const *x, unsigned int length, PDFType type,
					vector<double const *> const &par,
					vector<unsigned int> const &len,
					double const *lower, double const *upper) const
		{
			double t      = *x;
			double eta    = *par[0];
			double shape  = *par[1];

			return gamma_aft_log_density(t, eta, shape);
		}
		double SGAMMA::logDensity(double const *x, unsigned int length, PDFType type,
					vector<double const *> const &par,
					vector<unsigned int> const &len,
					double const *lower, double const *upper) const
		{
			double t      = *x;
			double eta    = *par[0];
			double shape  = *par[1];

			return gamma_aft_log_survival(t, eta, shape);
		}

		void DGAMMA::randomSample(double *x, unsigned int length,
					vector<double const *> const &par,
					vector<unsigned int> const &len,
					double const *lower, double const *upper,
					RNG *rng) const
		{
			// not implemented
		}
		void SGAMMA::randomSample(double *x, unsigned int length,
					vector<double const *> const &par,
					vector<unsigned int> const &len,
					double const *lower, double const *upper,
					RNG *rng) const
		{
			// not implemented
		}

		void DGAMMA::support(double *lower, double *upper, unsigned int length,
				vector<double const *> const &par,
				vector<unsigned int> const &len) const
		{
			for (unsigned int i = 0; i < length; ++i) {
				lower[i] = 0;
				upper[i] = JAGS_POSINF;
			}
		}
		void SGAMMA::support(double *lower, double *upper, unsigned int length,
				vector<double const *> const &par,
				vector<unsigned int> const &len) const
		{
			for (unsigned int i = 0; i < length; ++i) {
				lower[i] = 0;
				upper[i] = JAGS_POSINF;
			}
		}

		unsigned int DGAMMA::length(vector<unsigned int> const &len) const
		{
			return 1;
		}
		unsigned int SGAMMA::length(vector<unsigned int> const &len) const
		{
			return 1;
		}

		void DGAMMA::typicalValue(double *x, unsigned int length,
					vector<double const *> const &par,
					vector<unsigned int> const &len,
					double const *lower, double const *upper) const
		{
			// not implemented
		}
		void SGAMMA::typicalValue(double *x, unsigned int length,
					vector<double const *> const &par,
					vector<unsigned int> const &len,
					double const *lower, double const *upper) const
		{
			// not implemented
		}

		bool DGAMMA::isSupportFixed(vector<bool> const &fixmask) const
		{
			return true;
		}
		bool SGAMMA::isSupportFixed(vector<bool> const &fixmask) const
		{
			return true;
		}
	}
}

