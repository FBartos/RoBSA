#ifndef exp_H_
#define exp_H_

#include <function/ScalarFunction.h>

namespace jags {
  namespace surv {

    // effect sizes transformations
    class DEXP_fun : public ScalarFunction 
    {
      public:
        DEXP_fun();

        bool checkParameterValue(std::vector<double const *> const &args) const;
        double evaluate(std::vector<double const *> const &args) const;
    };

    class SEXP_fun : public ScalarFunction 
    {
      public:
        SEXP_fun();

        bool checkParameterValue(std::vector<double const *> const &args) const;
        double evaluate(std::vector<double const *> const &args) const;
    };
  }
}

#endif /* exp_H_ */
