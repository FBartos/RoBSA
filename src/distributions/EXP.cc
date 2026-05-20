
#include "EXP.h"
#include <util/nainf.h>
#include <cmath>
#include <rng/RNG.h>
#include <JRmath.h>
#include "../source/distributions.h"

using std::vector;

namespace jags {
	namespace surv {

		DEXP::DEXP() : VectorDist("exp_aft_event", 1) {}
		SEXP::SEXP() : VectorDist("exp_aft_cens_r", 1) {}

		bool DEXP::checkParameterLength(vector<unsigned long> const &len) const
		{
			return true;
		}
		bool SEXP::checkParameterLength(vector<unsigned long> const &len) const
		{
			return true;
		}

		bool DEXP::checkParameterValue(vector<double const *> const &par,
						vector<unsigned long> const &len) const
		{
			return true;
		}
		bool SEXP::checkParameterValue(vector<double const *> const &par,
						vector<unsigned long> const &len) const
		{
			return true;
		}

		double DEXP::logDensity(double const *x,  PDFType type,
					vector<double const *> const &par,
					vector<unsigned long> const &len) const
		{
			double t      = *x;
			double eta    = *par[0];

			return exp_aft_log_density(t, eta);
		}
		double SEXP::logDensity(double const *x,  PDFType type,
					vector<double const *> const &par,
					vector<unsigned long> const &len) const
		{
			double t      = *x;
			double eta    = *par[0];

			return exp_aft_log_survival(t, eta);
		}

		void DEXP::randomSample(double *x, 
					vector<double const *> const &par,
					vector<unsigned long> const &len,
					RNG *rng) const
		{
			// not implemented
		}
		void SEXP::randomSample(double *x, 
					vector<double const *> const &par,
					vector<unsigned long> const &len,
					RNG *rng) const
		{
			// not implemented
		}

		void DEXP::support(double *lower, double *upper, 
				vector<double const *> const &par,
				vector<unsigned long> const &len) const
		{
		    *lower = 0;
		    *upper = JAGS_POSINF;
		}
		void SEXP::support(double *lower, double *upper, 
				vector<double const *> const &par,
				vector<unsigned long> const &len) const
		{
		    *lower = 0;
		    *upper = JAGS_POSINF;
		}

		unsigned long DEXP::length(vector<unsigned long> const &len) const
		{
			return 1;
		}
		unsigned long SEXP::length(vector<unsigned long> const &len) const
		{
			return 1;
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

