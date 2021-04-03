#include "distribution_functions.h"
#include <util/nainf.h>
#include <JRmath.h>
#include <cmath>

using std::exp;
using std::sqrt;
using std::pow;
using std::log;

double exp_aft_log_density  (double t, double eta){
  return exp_aft_log_hazard(eta) + exp_aft_log_survival(t, eta);
}
double exp_aft_log_hazard   (double eta){
  return -eta;
}
double exp_aft_log_survival (double t, double eta){
  return -t * exp(-eta);
}

double weibull_aft_log_density  (double t, double eta, double shape){
  return weibull_aft_log_hazard(t, eta, shape) + weibull_aft_log_survival(t, eta, shape);
}
double weibull_aft_log_hazard   (double t, double eta, double shape){
  return log(shape) + (shape - 1) * log(t) - shape * eta;
}
double weibull_aft_log_survival (double t, double eta, double shape){
  return -pow(t, shape) * exp(-shape * eta);
}

double lnorm_aft_log_density  (double t, double eta, double sd){
  return dlnorm(t, eta, sd, true);
}
double lnorm_aft_log_hazard   (double t, double eta, double sd){
  return lnorm_aft_log_density(t, eta, sd) - lnorm_aft_log_survival(t, eta, sd);
}
double lnorm_aft_log_survival (double t, double eta, double sd){
  return plnorm(t, eta, sd, false, true);
}

double llogis_aft_log_density  (double t, double eta, double shape){
  return log(shape) - eta + (shape - 1) * (log(t) - eta) - 2 * log(1 + pow( t / exp(eta), shape));
}
double llogis_aft_log_hazard   (double t, double eta, double shape){
  return llogis_aft_log_density(t, eta, shape) - llogis_aft_log_survival(t, eta, shape);
}
double llogis_aft_log_survival (double t, double eta, double shape){
  return - log(1 + pow(t / exp(eta), shape));
}

double gamma_aft_log_density  (double t, double eta, double shape){
  return dgamma(t, shape, 1/ exp(-eta), true);
}
double gamma_aft_log_hazard   (double t, double eta, double shape){
  return gamma_aft_log_density(t, eta, shape) - gamma_aft_log_survival(t, eta, shape);
}
double gamma_aft_log_survival (double t, double eta, double shape){
  return pgamma(t, shape, 1/ exp(-eta), false, true);
}
