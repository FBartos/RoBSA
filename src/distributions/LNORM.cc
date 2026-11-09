#include "LNORM.h"
#include <util/nainf.h>
#include <cmath>
#include <rng/RNG.h>
#include <JRmath.h>
#include "../source/distributions.h"

using std::vector;

namespace jags {
	namespace surv {

		DLNORM::DLNORM() : VectorDist("lnorm_aft_event", 2) {}
		SLNORM::SLNORM() : VectorDist("lnorm_aft_cens_r", 2) {}

		bool DLNORM::checkParameterLength(vector<unsigned long> const &len) const
		{
			return true;
		}
		bool SLNORM::checkParameterLength(vector<unsigned long> const &len) const
		{
			return true;
		}

		bool DLNORM::checkParameterValue(vector<double const *> const &par,
						vector<unsigned long> const &len) const
		{
			return *par[1] > 0.0;
		}
		bool SLNORM::checkParameterValue(vector<double const *> const &par,
						vector<unsigned long> const &len) const
		{
			return *par[1] > 0.0;
		}

		double DLNORM::logDensity(double const *x,  PDFType type,
					vector<double const *> const &par,
					vector<unsigned long> const &len) const
		{
			double t      = *x;
			double eta    = *par[0];
			double sd     = *par[1];

			return lnorm_aft_log_density(t, eta, sd);
		}
		double SLNORM::logDensity(double const *x,  PDFType type,
					vector<double const *> const &par,
					vector<unsigned long> const &len) const
		{
			double t      = *x;
			double eta    = *par[0];
			double sd     = *par[1];

			return lnorm_aft_log_survival(t, eta, sd);
		}

		void DLNORM::randomSample(double *x, 
					vector<double const *> const &par,
					vector<unsigned long> const &len,
					RNG *rng) const
		{
			// not implemented
		}
		void SLNORM::randomSample(double *x, 
					vector<double const *> const &par,
					vector<unsigned long> const &len,
					RNG *rng) const
		{
			// not implemented
		}

		void DLNORM::support(double *lower, double *upper, 
				vector<double const *> const &par,
				vector<unsigned long> const &len) const
		{
		    *lower = 0;
		    *upper = JAGS_POSINF;
		}
		void SLNORM::support(double *lower, double *upper, 
				vector<double const *> const &par,
				vector<unsigned long> const &len) const
		{
		    *lower = 0;
		    *upper = JAGS_POSINF;
		}

		unsigned long DLNORM::length(vector<unsigned long> const &len) const
		{
			return 1;
		}
		unsigned long SLNORM::length(vector<unsigned long> const &len) const
		{
			return 1;
		}

		bool DLNORM::isSupportFixed(vector<bool> const &fixmask) const
		{
			return true;
		}
		bool SLNORM::isSupportFixed(vector<bool> const &fixmask) const
		{
			return true;
		}
	}
}

