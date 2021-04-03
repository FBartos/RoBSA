#ifndef distribution_functions_H_
#define distribution_functions_H_

double exp_aft_log_density      (double t, double eta);
double exp_aft_log_hazard       (double eta);
double exp_aft_log_survival     (double t, double eta);

double weibull_aft_log_density  (double t, double eta, double shape);
double weibull_aft_log_hazard   (double t, double eta, double shape);
double weibull_aft_log_survival (double t, double eta, double shape);

double lnorm_aft_log_density    (double t, double eta, double sd);
double lnorm_aft_log_hazard     (double t, double eta, double sd);
double lnorm_aft_log_survival   (double t, double eta, double sd);

double llogis_aft_log_density   (double t, double eta, double shape);
double llogis_aft_log_hazard    (double t, double eta, double shape);
double llogis_aft_log_survival  (double t, double eta, double shape);

double gamma_aft_log_density    (double t, double eta, double shape);
double gamma_aft_log_hazard     (double t, double eta, double shape);
double gamma_aft_log_survival   (double t, double eta, double shape);
#endif /* distribution_functions_H_ */
