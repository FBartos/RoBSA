
#include "LLOGIS.h"
#include <util/nainf.h>
#include <cmath>
#include <rng/RNG.h>
#include <JRmath.h>
#include "../source/distributions.h"

using std::vector;

namespace jags {
	namespace surv {

		DLLOGIS::DLLOGIS() : VectorDist("llogis_aft_event", 2) {}
		SLLOGIS::SLLOGIS() : VectorDist("llogis_aft_cens_r", 2) {}

		bool DLLOGIS::checkParameterLength(vector<unsigned int> const &len) const
		{
			return true;
		}
		bool SLLOGIS::checkParameterLength(vector<unsigned int> const &len) const
		{
			return true;
		}

		bool DLLOGIS::checkParameterValue(vector<double const *> const &par,
						vector<unsigned int> const &len) const
		{
			return *par[1] > 0.0;
		}
		bool SLLOGIS::checkParameterValue(vector<double const *> const &par,
						vector<unsigned int> const &len) const
		{
			return *par[1] > 0.0;
		}

		double DLLOGIS::logDensity(double const *x, unsigned int length, PDFType type,
					vector<double const *> const &par,
					vector<unsigned int> const &len,
					double const *lower, double const *upper) const
		{
			double t     = *x;
			double eta   = *par[0];
			double shape = *par[1];
			return llogis_aft_log_density(t, eta, shape);
		}
		double SLLOGIS::logDensity(double const *x, unsigned int length, PDFType type,
					vector<double const *> const &par,
					vector<unsigned int> const &len,
					double const *lower, double const *upper) const
		{
			double t     = *x;
			double eta   = *par[0];
			double shape = *par[1];
			return llogis_aft_log_survival(t, eta, shape);
		}

		void DLLOGIS::randomSample(double *x, unsigned int length,
					vector<double const *> const &par,
					vector<unsigned int> const &len,
					double const *lower, double const *upper,
					RNG *rng) const
		{
			// not implemented
		}
		void SLLOGIS::randomSample(double *x, unsigned int length,
					vector<double const *> const &par,
					vector<unsigned int> const &len,
					double const *lower, double const *upper,
					RNG *rng) const
		{
			// not implemented
		}

		void DLLOGIS::support(double *lower, double *upper, unsigned int length,
				vector<double const *> const &par,
				vector<unsigned int> const &len) const
		{
			for (unsigned int i = 0; i < length; ++i) {
				lower[i] = 0;
				upper[i] = JAGS_POSINF;
			}
		}
		void SLLOGIS::support(double *lower, double *upper, unsigned int length,
				vector<double const *> const &par,
				vector<unsigned int> const &len) const
		{
			for (unsigned int i = 0; i < length; ++i) {
				lower[i] = 0;
				upper[i] = JAGS_POSINF;
			}
		}

		unsigned int DLLOGIS::length(vector<unsigned int> const &len) const
		{
			return 1;
		}
		unsigned int SLLOGIS::length(vector<unsigned int> const &len) const
		{
			return 1;
		}

		void DLLOGIS::typicalValue(double *x, unsigned int length,
					vector<double const *> const &par,
					vector<unsigned int> const &len,
					double const *lower, double const *upper) const
		{
			// not implemented
		}
		void SLLOGIS::typicalValue(double *x, unsigned int length,
					vector<double const *> const &par,
					vector<unsigned int> const &len,
					double const *lower, double const *upper) const
		{
			// not implemented
		}

		bool DLLOGIS::isSupportFixed(vector<bool> const &fixmask) const
		{
			return true;
		}
		bool SLLOGIS::isSupportFixed(vector<bool> const &fixmask) const
		{
			return true;
		}
	}
}
