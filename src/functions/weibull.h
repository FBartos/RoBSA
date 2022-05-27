#ifndef weibull_H_
#define weibull_H_

#include <function/ScalarFunction.h>

namespace jags {
  namespace surv {

    // effect sizes transformations
    class DWEIBULL_fun : public ScalarFunction 
    {
      public:
        DWEIBULL_fun();

        bool checkParameterValue(std::vector<double const *> const &args) const;
        double evaluate(std::vector<double const *> const &args) const;
    };

    class SWEIBULL_fun : public ScalarFunction 
    {
      public:
        SWEIBULL_fun();

        bool checkParameterValue(std::vector<double const *> const &args) const;
        double evaluate(std::vector<double const *> const &args) const;
    };
  }
}

#endif /* weibull_H_ */
