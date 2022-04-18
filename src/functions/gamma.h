#ifndef gamma_H_
#define gamma_H_

#include <function/ScalarFunction.h>

namespace jags {
  namespace surv {

    // effect sizes transformations
    class DGAMMA_fun : public ScalarFunction 
    {
      public:
        DGAMMA_fun();

        bool checkParameterValue(std::vector<double const *> const &args) const;
        double evaluate(std::vector<double const *> const &args) const;
    };

    class SGAMMA_fun : public ScalarFunction 
    {
      public:
        SGAMMA_fun();

        bool checkParameterValue(std::vector<double const *> const &args) const;
        double evaluate(std::vector<double const *> const &args) const;
    };
  }
}

#endif /* gamma_H_ */
