#ifndef INTEGRAL_H
#define INTEGRAL_H

#include<cmath>
#include"elements.h"

double integral_q(double (*f)(double ,double, double, double, double, double, double), double m0, double q_max, double z, double lambda_e, double lambda_Pi, double gamma_1, double gamma_2);

double integral_z_for_e(double lambda_e, double lambda_Pi, double gamma_1, double gamma_2);

double integral_z_for_P(double lambda_e, double lambda_Pi, double gamma_1, double gamma_2);

double integral_z_for_pi1(double lambda_e, double lambda_Pi, double gamma_1, double gamma_2);

double integral_z_for_pi2(double lambda_e, double lambda_Pi, double gamma_1, double gamma_2);



#endif