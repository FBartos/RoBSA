
#include "LNORM.h"
#include <util/nainf.h>
#include <cmath>
#include <rng/RNG.h>
#include <JRmath.h>
#include "distribution_functions.h"

using std::vector;

namespace jags {
	namespace surv {

		DLNORM::DLNORM() : VectorDist("lnorm_aft_event", 2) {}
		SLNORM::SLNORM() : VectorDist("lnorm_aft_rcent", 2) {}

		bool DLNORM::checkParameterLength(vector<unsigned int> const &len) const
		{
			return true;
		}
		bool SLNORM::checkParameterLength(vector<unsigned int> const &len) const
		{
			return true;
		}

		bool DLNORM::checkParameterValue(vector<double const *> const &par,
						vector<unsigned int> const &len) const
		{
			return *par[1] > 0.0;
		}
		bool SLNORM::checkParameterValue(vector<double const *> const &par,
						vector<unsigned int> const &len) const
		{
			return *par[1] > 0.0;
		}

		double DLNORM::logDensity(double const *x, unsigned int length, PDFType type,
					vector<double const *> const &par,
					vector<unsigned int> const &len,
					double const *lower, double const *upper) const
		{
			double t      = *x;
			double eta    = *par[0];
			double sd     = *par[1];

			return lnorm_aft_log_density(t, eta, sd);
		}
		double SLNORM::logDensity(double const *x, unsigned int length, PDFType type,
					vector<double const *> const &par,
					vector<unsigned int> const &len,
					double const *lower, double const *upper) const
		{
			double t      = *x;
			double eta    = *par[0];
			double sd     = *par[1];

			return lnorm_aft_log_survival(t, eta, sd);
		}

		void DLNORM::randomSample(double *x, unsigned int length,
					vector<double const *> const &par,
					vector<unsigned int> const &len,
					double const *lower, double const *upper,
					RNG *rng) const
		{
			// not implemented
		}
		void SLNORM::randomSample(double *x, unsigned int length,
					vector<double const *> const &par,
					vector<unsigned int> const &len,
					double const *lower, double const *upper,
					RNG *rng) const
		{
			// not implemented
		}

		void DLNORM::support(double *lower, double *upper, unsigned int length,
				vector<double const *> const &par,
				vector<unsigned int> const &len) const
		{
			for (unsigned int i = 0; i < length; ++i) {
				lower[i] = 0;
				upper[i] = JAGS_POSINF;
			}
		}
		void SLNORM::support(double *lower, double *upper, unsigned int length,
				vector<double const *> const &par,
				vector<unsigned int> const &len) const
		{
			for (unsigned int i = 0; i < length; ++i) {
				lower[i] = 0;
				upper[i] = JAGS_POSINF;
			}
		}

		unsigned int DLNORM::length(vector<unsigned int> const &len) const
		{
			return 1;
		}
		unsigned int SLNORM::length(vector<unsigned int> const &len) const
		{
			return 1;
		}

		void DLNORM::typicalValue(double *x, unsigned int length,
					vector<double const *> const &par,
					vector<unsigned int> const &len,
					double const *lower, double const *upper) const
		{
			// not implemented
		}
		void SLNORM::typicalValue(double *x, unsigned int length,
					vector<double const *> const &par,
					vector<unsigned int> const &len,
					double const *lower, double const *upper) const
		{
			// not implemented
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

