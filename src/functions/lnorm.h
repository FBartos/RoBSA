#ifndef lnorm_H_
#define lnorm_H_

#include <function/ScalarFunction.h>

namespace jags {
  namespace surv {

    // effect sizes transformations
    class DLNORM_fun : public ScalarFunction 
    {
      public:
        DLNORM_fun();

        bool checkParameterValue(std::vector<double const *> const &args) const;
        double evaluate(std::vector<double const *> const &args) const;
    };

    class SLNORM_fun : public ScalarFunction 
    {
      public:
        SLNORM_fun();

        bool checkParameterValue(std::vector<double const *> const &args) const;
        double evaluate(std::vector<double const *> const &args) const;
    };
  }
}

#endif /* lnorm_H_ */
