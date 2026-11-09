#ifndef LLOGIS_H_
#define LLOGIS_H_
#include <distribution/VectorDist.h>

namespace jags {
  namespace surv { // module namespace

    class DLLOGIS : public VectorDist
    {
      public:
        DLLOGIS();
        
      double logDensity(double const *x,  PDFType tpye, 
            std::vector<double const *> const &parameters,
            std::vector<unsigned long> const &lengths) const override;
      void randomSample(double *x, 
            std::vector<double const *> const &parameters,
            std::vector<unsigned long> const &lengths, RNG *rng) const override;
      bool checkParameterValue(std::vector<double const *> const &parameters,
                              std::vector<unsigned long> const &lengths) const override;
      bool checkParameterLength(std::vector<unsigned long> const &lengths) const override;
      unsigned long length(std::vector<unsigned long> const &dim) const override;
      void support(double *lower, double *upper, 
            std::vector<double const *> const &parameters,
            std::vector<unsigned long> const &lengths) const override;
      bool isSupportFixed(std::vector<bool> const &fixmask) const override;
    };

    class SLLOGIS : public VectorDist
    {
      public:
        SLLOGIS();
        
      double logDensity(double const *x,  PDFType tpye, 
            std::vector<double const *> const &parameters,
            std::vector<unsigned long> const &lengths) const override;
      void randomSample(double *x, 
            std::vector<double const *> const &parameters,
            std::vector<unsigned long> const &lengths, RNG *rng) const override;
      bool checkParameterValue(std::vector<double const *> const &parameters,
                              std::vector<unsigned long> const &lengths) const override;
      bool checkParameterLength(std::vector<unsigned long> const &lengths) const override;
      unsigned long length(std::vector<unsigned long> const &dim) const override;
      void support(double *lower, double *upper, 
            std::vector<double const *> const &parameters,
            std::vector<unsigned long> const &lengths) const override;
      bool isSupportFixed(std::vector<bool> const &fixmask) const override;
    };

  }
}

#endif /* LLOGIS_H_ */
