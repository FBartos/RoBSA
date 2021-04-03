
#include "LLOGIST.h"
#include <util/nainf.h>
#include <cmath>
#include <rng/RNG.h>
#include <JRmath.h>
#include "distribution_functions.h"

using std::vector;

namespace jags {
	namespace surv {

		DLLOGIST::DLLOGIST() : VectorDist("llogis_aft_event", 2) {}
		SLLOGIST::SLLOGIST() : VectorDist("llogis_aft_rcent", 2) {}

		bool DLLOGIST::checkParameterLength(vector<unsigned int> const &len) const
		{
			return true;
		}
		bool SLLOGIST::checkParameterLength(vector<unsigned int> const &len) const
		{
			return true;
		}

		bool DLLOGIST::checkParameterValue(vector<double const *> const &par,
						vector<unsigned int> const &len) const
		{
			return *par[1] > 0.0;
		}
		bool SLLOGIST::checkParameterValue(vector<double const *> const &par,
						vector<unsigned int> const &len) const
		{
			return *par[1] > 0.0;
		}

		double DLLOGIST::logDensity(double const *x, unsigned int length, PDFType type,
					vector<double const *> const &par,
					vector<unsigned int> const &len,
					double const *lower, double const *upper) const
		{
			double t     = *x;
			double eta   = *par[0];
			double shape = *par[1];
			return llogis_aft_log_density(t, eta, shape);
		}
		double SLLOGIST::logDensity(double const *x, unsigned int length, PDFType type,
					vector<double const *> const &par,
					vector<unsigned int> const &len,
					double const *lower, double const *upper) const
		{
			double t     = *x;
			double eta   = *par[0];
			double shape = *par[1];
			return llogis_aft_log_survival(t, eta, shape);
		}

		void DLLOGIST::randomSample(double *x, unsigned int length,
					vector<double const *> const &par,
					vector<unsigned int> const &len,
					double const *lower, double const *upper,
					RNG *rng) const
		{
			// not implemented
		}
		void SLLOGIST::randomSample(double *x, unsigned int length,
					vector<double const *> const &par,
					vector<unsigned int> const &len,
					double const *lower, double const *upper,
					RNG *rng) const
		{
			// not implemented
		}

		void DLLOGIST::support(double *lower, double *upper, unsigned int length,
				vector<double const *> const &par,
				vector<unsigned int> const &len) const
		{
			for (unsigned int i = 0; i < length; ++i) {
				lower[i] = 0;
				upper[i] = JAGS_POSINF;
			}
		}
		void SLLOGIST::support(double *lower, double *upper, unsigned int length,
				vector<double const *> const &par,
				vector<unsigned int> const &len) const
		{
			for (unsigned int i = 0; i < length; ++i) {
				lower[i] = 0;
				upper[i] = JAGS_POSINF;
			}
		}

		unsigned int DLLOGIST::length(vector<unsigned int> const &len) const
		{
			return 1;
		}
		unsigned int SLLOGIST::length(vector<unsigned int> const &len) const
		{
			return 1;
		}

		void DLLOGIST::typicalValue(double *x, unsigned int length,
					vector<double const *> const &par,
					vector<unsigned int> const &len,
					double const *lower, double const *upper) const
		{
			// not implemented
		}
		void SLLOGIST::typicalValue(double *x, unsigned int length,
					vector<double const *> const &par,
					vector<unsigned int> const &len,
					double const *lower, double const *upper) const
		{
			// not implemented
		}

		bool DLLOGIST::isSupportFixed(vector<bool> const &fixmask) const
		{
			return true;
		}
		bool SLLOGIST::isSupportFixed(vector<bool> const &fixmask) const
		{
			return true;
		}
	}
}
