#ifndef ELEMENTS_H
#define ELEMENTS_H
#include<iostream>
#include<cmath>

double dedqdz_withm0(double m0, double z, double lambda_e, double lambda_Pi, double gamma_1, double gamma_2, double q);

double dPdqdz_withm0(double m0, double z, double lambda_e, double lambda_Pi, double gamma_1, double gamma_2, double q);

double dpi1dqdz_withm0(double m0, double z, double lambda_e, double lambda_Pi, double gamma_1, double gamma_2, double q);

double dpi2dqdz_withm0(double m0, double z, double lambda_e, double lambda_Pi, double gamma_1, double gamma_2, double q);

double analytic_dedz_withoutm0(double z, double lambda_e, double lambda_Pi, double gamma_1, double gamma_2);

double analytic_dPdz_withoutm0(double z, double lambda_e, double lambda_Pi, double gamma_1, double gamma_2);

double analytic_dpi1dz_withoutm0(double z, double lambda_e, double lambda_Pi, double gamma_1, double gamma_2);

double analytic_dpi2dz_withoutm0(double z, double lambda_e, double lambda_Pi, double gamma_1, double gamma_2);


#endif