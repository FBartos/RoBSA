
#include "EXP.h"
#include <util/nainf.h>
#include <cmath>
#include <rng/RNG.h>
#include <JRmath.h>
#include "distribution_functions.h"

using std::vector;

namespace jags {
	namespace surv {

		DEXP::DEXP() : VectorDist("exp_aft_event", 1) {}
		SEXP::SEXP() : VectorDist("exp_aft_cens_r", 1) {}

		bool DEXP::checkParameterLength(vector<unsigned int> const &len) const
		{
			return true;
		}
		bool SEXP::checkParameterLength(vector<unsigned int> const &len) const
		{
			return true;
		}

		bool DEXP::checkParameterValue(vector<double const *> const &par,
						vector<unsigned int> const &len) const
		{
			return true;
		}
		bool SEXP::checkParameterValue(vector<double const *> const &par,
						vector<unsigned int> const &len) const
		{
			return true;
		}

		double DEXP::logDensity(double const *x, unsigned int length, PDFType type,
					vector<double const *> const &par,
					vector<unsigned int> const &len,
					double const *lower, double const *upper) const
		{
			double t      = *x;
			double eta    = *par[0];

			return exp_aft_log_density(t, eta);
		}
		double SEXP::logDensity(double const *x, unsigned int length, PDFType type,
					vector<double const *> const &par,
					vector<unsigned int> const &len,
					double const *lower, double const *upper) const
		{
			double t      = *x;
			double eta    = *par[0];

			return exp_aft_log_survival(t, eta);
		}

		void DEXP::randomSample(double *x, unsigned int length,
					vector<double const *> const &par,
					vector<unsigned int> const &len,
					double const *lower, double const *upper,
					RNG *rng) const
		{
			// not implemented
		}
		void SEXP::randomSample(double *x, unsigned int length,
					vector<double const *> const &par,
					vector<unsigned int> const &len,
					double const *lower, double const *upper,
					RNG *rng) const
		{
			// not implemented
		}

		void DEXP::support(double *lower, double *upper, unsigned int length,
				vector<double const *> const &par,
				vector<unsigned int> const &len) const
		{
			for (unsigned int i = 0; i < length; ++i) {
				lower[i] = 0;
				upper[i] = JAGS_POSINF;
			}
		}
		void SEXP::support(double *lower, double *upper, unsigned int length,
				vector<double const *> const &par,
				vector<unsigned int> const &len) const
		{
			for (unsigned int i = 0; i < length; ++i) {
				lower[i] = 0;
				upper[i] = JAGS_POSINF;
			}
		}

		unsigned int DEXP::length(vector<unsigned int> const &len) const
		{
			return 1;
		}
		unsigned int SEXP::length(vector<unsigned int> const &len) const
		{
			return 1;
		}


		void DEXP::typicalValue(double *x, unsigned int length,
					vector<double const *> const &par,
					vector<unsigned int> const &len,
					double const *lower, double const *upper) const
		{
			// not implemented
		}
		void SEXP::typicalValue(double *x, unsigned int length,
					vector<double const *> const &par,
					vector<unsigned int> const &len,
					double const *lower, double const *upper) const
		{
			// not implemented
		}

		bool DEXP::isSupportFixed(vector<bool> const &fixmask) const
		{
			return true;
		}
		bool SEXP::isSupportFixed(vector<bool> const &fixmask) const
		{
			return true;
		}
	}
}

