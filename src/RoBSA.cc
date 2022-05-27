#include <module/Module.h>
#include "distributions/EXP.h"
#include "distributions/WEIBULL.h"
#include "distributions/LNORM.h"
#include "distributions/LLOGIS.h"
#include "distributions/GAMMA.h"

#include "functions/exp.h"
#include "functions/weibull.h"
#include "functions/lnorm.h"
#include "functions/llogis.h"
#include "functions/gamma.h"

namespace jags {
  namespace surv { // module namespace

    // JAGS module class
    class SurvModule : public Module {
      public:
        SurvModule();
        ~SurvModule();
    };

    // constructor (executed when loading the module)
    SurvModule::SurvModule() : Module("RoBSA"){

      // density functions for observed events
      insert(new DEXP);
      insert(new DWEIBULL);
      insert(new DLNORM);
      insert(new DLLOGIS);
      insert(new DGAMMA);

      // density functions for right truncated events (survival functions)
      insert(new SEXP);
      insert(new SWEIBULL);
      insert(new SLNORM);
      insert(new SLLOGIS);
      insert(new SGAMMA);

      // likelihood functions
      insert(new DEXP_fun);
      insert(new DWEIBULL_fun);
      insert(new DLNORM_fun);
      insert(new DLLOGIS_fun);
      insert(new DGAMMA_fun);

      insert(new SEXP_fun);
      insert(new SWEIBULL_fun);
      insert(new SLNORM_fun);
      insert(new SLLOGIS_fun);
      insert(new SGAMMA_fun);
    }

    // destructor (executed when unloading the module)
    SurvModule::~SurvModule() {
      std::vector<Function*> const &fvec = functions();
      for (unsigned int i = 0; i < fvec.size(); ++i) {
        delete fvec[i];
      }
      std::vector<Distribution*> const &dvec = distributions();
      for (unsigned int i=0;i<dvec.size();++i) {
        delete dvec[i];
      }
    }

  }
}

jags::surv::SurvModule _surv_module;

