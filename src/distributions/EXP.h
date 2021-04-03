#ifndef EXP_H_
#define EXP_H_
#include <distribution/VectorDist.h>

namespace jags {
  namespace surv { // module namespace

    class DEXP : public VectorDist
    {
      public:
        DEXP();
        
      double logDensity(double const *x, unsigned int length, PDFType tpye, 
            std::vector<double const *> const &parameters,
            std::vector<unsigned int> const &lengths,
            double const *lower, double const *upper) const;
      void randomSample(double *x, unsigned int length,
            std::vector<double const *> const &parameters,
            std::vector<unsigned int> const &lengths,
            double const *lower, double const *upper, RNG *rng) const;
      void typicalValue(double *x, unsigned int length,
            std::vector<double const *> const &par,
            std::vector<unsigned int> const &lengths,
            double const *lower, double const *upper) const;
      bool checkParameterValue(std::vector<double const *> const &parameters,
                              std::vector<unsigned int> const &lengths) const;
      bool checkParameterLength(std::vector<unsigned int> const &lengths) const;
      unsigned int length(std::vector<unsigned int> const &dim) const;
      void support(double *lower, double *upper, unsigned int length,
            std::vector<double const *> const &parameters,
            std::vector<unsigned int> const &lengths) const;
      bool isSupportFixed(std::vector<bool> const &fixmask) const;
    };

    class SEXP : public VectorDist
    {
      public:
        SEXP();
        
      double logDensity(double const *x, unsigned int length, PDFType tpye, 
            std::vector<double const *> const &parameters,
            std::vector<unsigned int> const &lengths,
            double const *lower, double const *upper) const;
      void randomSample(double *x, unsigned int length,
            std::vector<double const *> const &parameters,
            std::vector<unsigned int> const &lengths,
            double const *lower, double const *upper, RNG *rng) const;
      void typicalValue(double *x, unsigned int length,
            std::vector<double const *> const &par,
            std::vector<unsigned int> const &lengths,
            double const *lower, double const *upper) const;
      bool checkParameterValue(std::vector<double const *> const &parameters,
                              std::vector<unsigned int> const &lengths) const;
      bool checkParameterLength(std::vector<unsigned int> const &lengths) const;
      unsigned int length(std::vector<unsigned int> const &dim) const;
      void support(double *lower, double *upper, unsigned int length,
            std::vector<double const *> const &parameters,
            std::vector<unsigned int> const &lengths) const;
      bool isSupportFixed(std::vector<bool> const &fixmask) const;
    };

  }
}

#endif /* EXP_H_ */