
#include "GAMMA.h"
#include <util/nainf.h>
#include <cmath>
#include <rng/RNG.h>
#include <JRmath.h>
#include "../source/distributions.h"

using std::vector;

namespace jags {
	namespace surv {

		DGAMMA::DGAMMA() : VectorDist("gamma_aft_event", 2) {}
		SGAMMA::SGAMMA() : VectorDist("gamma_aft_cens_r", 2) {}

		bool DGAMMA::checkParameterLength(vector<unsigned long> const &len) const
		{
			return true;
		}
		bool SGAMMA::checkParameterLength(vector<unsigned long> const &len) const
		{
			return true;
		}

		bool DGAMMA::checkParameterValue(vector<double const *> const &par,
						vector<unsigned long> const &len) const
		{
			return *par[1] > 0.0;
		}
		bool SGAMMA::checkParameterValue(vector<double const *> const &par,
						vector<unsigned long> const &len) const
		{
			return *par[1] > 0.0;
		}

		double DGAMMA::logDensity(double const *x,  PDFType type,
					vector<double const *> const &par,
					vector<unsigned long> const &len) const
		{
			double t      = *x;
			double eta    = *par[0];
			double shape  = *par[1];

			return gamma_aft_log_density(t, eta, shape);
		}
		double SGAMMA::logDensity(double const *x,  PDFType type,
					vector<double const *> const &par,
					vector<unsigned long> const &len) const
		{
			double t      = *x;
			double eta    = *par[0];
			double shape  = *par[1];

			return gamma_aft_log_survival(t, eta, shape);
		}

		void DGAMMA::randomSample(double *x, 
					vector<double const *> const &par,
					vector<unsigned long> const &len,
					RNG *rng) const
		{
			// not implemented
		}
		void SGAMMA::randomSample(double *x, 
					vector<double const *> const &par,
					vector<unsigned long> const &len,
					RNG *rng) const
		{
			// not implemented
		}

		void DGAMMA::support(double *lower, double *upper, 
				vector<double const *> const &par,
				vector<unsigned long> const &len) const
		{
		    *lower = 0;
		    *upper = JAGS_POSINF;
		}
		void SGAMMA::support(double *lower, double *upper, 
				vector<double const *> const &par,
				vector<unsigned long> const &len) const
		{
		    *lower = 0;
		    *upper = JAGS_POSINF;
		}

		unsigned long DGAMMA::length(vector<unsigned long> const &len) const
		{
			return 1;
		}
		unsigned long SGAMMA::length(vector<unsigned long> const &len) const
		{
			return 1;
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

