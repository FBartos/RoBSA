#ifndef llogis_H_
#define llogis_H_

#include <function/ScalarFunction.h>

namespace jags {
  namespace surv {

    // effect sizes transformations
    class DLLOGIS_fun : public ScalarFunction 
    {
      public:
        DLLOGIS_fun();

        bool checkParameterValue(std::vector<double const *> const &args) const;
        double evaluate(std::vector<double const *> const &args) const;
    };

    class SLLOGIS_fun : public ScalarFunction 
    {
      public:
        SLLOGIS_fun();

        bool checkParameterValue(std::vector<double const *> const &args) const;
        double evaluate(std::vector<double const *> const &args) const;
    };
  }
}

#endif /* llogis_H_ */
