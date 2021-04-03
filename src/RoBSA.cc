#include <module/Module.h>
#include "distributions/EXP.h"
#include "distributions/WEIBULL.h"
#include "distributions/LNORM.h"
#include "distributions/LLOGIST.h"
#include "distributions/GAMMA.h"

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
      insert(new DLLOGIST);
      insert(new DGAMMA);

      // density functions for right truncated events (survival functions)
      insert(new SEXP);
      insert(new SWEIBULL);
      insert(new SLNORM);
      insert(new SLLOGIST);
      insert(new SGAMMA);


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

