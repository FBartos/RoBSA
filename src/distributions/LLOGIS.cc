
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

		bool DLLOGIS::checkParameterLength(vector<unsigned long> const &len) const
		{
			return true;
		}
		bool SLLOGIS::checkParameterLength(vector<unsigned long> const &len) const
		{
			return true;
		}

		bool DLLOGIS::checkParameterValue(vector<double const *> const &par,
						vector<unsigned long> const &len) const
		{
			return *par[1] > 0.0;
		}
		bool SLLOGIS::checkParameterValue(vector<double const *> const &par,
						vector<unsigned long> const &len) const
		{
			return *par[1] > 0.0;
		}

		double DLLOGIS::logDensity(double const *x,  PDFType type,
					vector<double const *> const &par,
					vector<unsigned long> const &len) const
		{
			double t     = *x;
			double eta   = *par[0];
			double shape = *par[1];
			return llogis_aft_log_density(t, eta, shape);
		}
		double SLLOGIS::logDensity(double const *x,  PDFType type,
					vector<double const *> const &par,
					vector<unsigned long> const &len) const
		{
			double t     = *x;
			double eta   = *par[0];
			double shape = *par[1];
			return llogis_aft_log_survival(t, eta, shape);
		}

		void DLLOGIS::randomSample(double *x, 
					vector<double const *> const &par,
					vector<unsigned long> const &len,
					RNG *rng) const
		{
			// not implemented
		}
		void SLLOGIS::randomSample(double *x, 
					vector<double const *> const &par,
					vector<unsigned long> const &len,
					RNG *rng) const
		{
			// not implemented
		}

		void DLLOGIS::support(double *lower, double *upper, 
				vector<double const *> const &par,
				vector<unsigned long> const &len) const
		{
		    *lower = 0;
		    *upper = JAGS_POSINF;
		}
		void SLLOGIS::support(double *lower, double *upper, 
				vector<double const *> const &par,
				vector<unsigned long> const &len) const
		{
		    *lower = 0;
		    *upper = JAGS_POSINF;
		}

		unsigned long DLLOGIS::length(vector<unsigned long> const &len) const
		{
			return 1;
		}
		unsigned long SLLOGIS::length(vector<unsigned long> const &len) const
		{
			return 1;
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
