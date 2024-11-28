#include<iostream>
#include<cmath>
#include"elements.h"

const double PI = std::acos(-1);

double dedqdz_withm0(double m0, double z, double lambda_e, double lambda_Pi, double gamma_1, double gamma_2, double q)
{
    if (q < 1e-8){q = 1e-8;}
    double q_sq = q*q;
    double q0 = std::sqrt(q_sq+m0*m0);
    double inside_exp = -lambda_Pi*q_sq/q0 - lambda_e*q0 - (gamma_1+gamma_2)*(1.-3.*z*z)*q_sq/(2.*q0);
    double exponential_term = std::exp(inside_exp);
    double inside_bessel0 = (gamma_1-gamma_2)*(z*z-1.)*q_sq/(2.*q0);
    //std::cout << inside_bessel0 << std::endl;
    try 
        {
        double bessel0_term =  std::cyl_bessel_i(0, abs(inside_bessel0));
        } catch (const std::domain_error& e) {
            std::cout << "error of nan in bessel" ;
            exit(-1);
                                                }
    double bessel0_term =  std::cyl_bessel_i(0, abs(inside_bessel0));
    double res = exponential_term*bessel0_term*q0*q_sq/(2.*PI*PI);
    return res;
    
}

double dPdqdz_withm0(double m0, double z, double lambda_e, double lambda_Pi, double gamma_1, double gamma_2, double q)
{
    if (q < 1e-8){q = 1e-8;}
    double q_sq = q*q;
    double q0 = std::sqrt(q_sq+m0*m0);
    double inside_exp = -lambda_Pi*q_sq/q0 - lambda_e*q0 - (gamma_1+gamma_2)*(1.-3.*z*z)*q_sq/(2.*q0);
    double exponential_term = std::exp(inside_exp);
    double inside_bessel0 = (gamma_1-gamma_2)*(z*z-1.)*q_sq/(2.*q0);
    double bessel0_term =  std::cyl_bessel_i(0, std::abs(inside_bessel0));
    double res = exponential_term*bessel0_term*q_sq*q_sq/(q0*2.*PI*PI);
    return res/3.;
}

double dpi1dqdz_withm0(double m0, double z, double lambda_e, double lambda_Pi, double gamma_1, double gamma_2, double q)
{
    if (q < 1e-8){q = 1e-8;}
    double q_sq = q*q;
    double q0 = std::sqrt(q_sq+m0*m0);
    double inside_exp = -lambda_Pi*q_sq/q0 - lambda_e*q0 - (gamma_1+gamma_2)*(1.-3.*z*z)*q_sq/(2.*q0);
    double exponential_term = std::exp(inside_exp);
    double inside_bessel0 = (gamma_1-gamma_2)*(z*z-1.)*q_sq/(2.*q0);
    double bessel0_term =  std::cyl_bessel_i(0, std::abs(inside_bessel0));
    double inside_bessel1 = (gamma_1-gamma_2)*q_sq*(1.-z*z)/(2.*q0);
    double sign_term = 1.;
    if(inside_bessel1<0){sign_term = -1.;}
    double bessel1_term =  sign_term*std::cyl_bessel_i(1, std::abs(inside_bessel1));
    double res = exponential_term*q_sq*q_sq*((-1.+3.*z*z)*bessel0_term - 3.*(-1.+z*z)*bessel1_term)/(6.*q0*2.*PI*PI);
    return -res;
}

double dpi2dqdz_withm0(double m0, double z, double lambda_e, double lambda_Pi, double gamma_1, double gamma_2, double q)
{
    if (q < 1e-8){q = 1e-8;}
    double q_sq = q*q;
    double q0 = std::sqrt(q_sq+m0*m0);
    double inside_exp = -lambda_Pi*q_sq/q0 - lambda_e*q0 - (gamma_1+gamma_2)*(1.-3.*z*z)*q_sq/(2.*q0);
    double exponential_term = std::exp(inside_exp);
    double inside_bessel0 = (gamma_1-gamma_2)*(z*z-1.)*q_sq/(2.*q0);
    double bessel0_term =  std::cyl_bessel_i(0, std::abs(inside_bessel0));
    double inside_bessel1 = (gamma_1-gamma_2)*q_sq*(1.-z*z)/(2.*q0);
    double sign_term = 1.;
    if(inside_bessel1<0){sign_term = -1.;}
    double bessel1_term =  sign_term*std::cyl_bessel_i(1, std::abs(inside_bessel1));
    double res = exponential_term*q_sq*q_sq*((-1.+3.*z*z)*bessel0_term + 3.*(-1.+z*z)*bessel1_term)/(6.*q0*2.*PI*PI);
    return -res;
}

double analytic_dedz_withoutm0(double z, double lambda_e, double lambda_Pi, double gamma_1, double gamma_2)
{
    double term_num_1 = -2.*(lambda_e + lambda_Pi) + gamma_1*(-1.+3.*z*z)+ gamma_2*(-1.+3.*z*z);
    double term_num_2 = 8.*(lambda_e + lambda_Pi)*(lambda_e + lambda_Pi) 
                        -8.*gamma_1*(lambda_e + lambda_Pi)*(-1.+3.*z*z)
                        -8.*gamma_2*(lambda_e + lambda_Pi)*(-1.+3.*z*z)
                        +2.*gamma_1*gamma_2*(-1.-6.*z*z+15.*z*z*z*z)
                        +gamma_1*gamma_1*(5.-18.*z*z+21.*z*z*z*z)
                        +gamma_2*gamma_2*(5.-18.*z*z+21.*z*z*z*z);
    double term_num = 24.*term_num_1*term_num_2;
    double term_denom = PI*PI*std::pow((-1.*(gamma_1-gamma_2)*(gamma_1-gamma_2)*(-1.+z*z)*(-1.+z*z)
                        + std::pow((gamma_1+gamma_2+2.*(lambda_e+lambda_Pi)-3.*gamma_1*z*z-3.*gamma_2*z*z),2.)),7./2.);
    //std::cout <<term_num<<std::endl;
    double res = -term_num/term_denom;
    //std::cout << lambda_e << " " << lambda_Pi << " " 
    //<< gamma_1 << " " << gamma_2 << " " << z << " "
    //<<term_num<< " " 
    //<< term_denom  << "  " << res <<std::endl;
    return res;
}

double analytic_dPdz_withoutm0(double z, double lambda_e, double lambda_Pi, double gamma_1, double gamma_2)
{
    double term_num_1 = -2*(lambda_e + lambda_Pi) + gamma_1*(-1+3*z*z)+ gamma_2*(-1+3*z*z);
    double term_num_2 = 8*(lambda_e + lambda_Pi)*(lambda_e + lambda_Pi) 
                        -8*gamma_1*(lambda_e + lambda_Pi)*(-1+3*z*z)
                        -8*gamma_2*(lambda_e + lambda_Pi)*(-1+3*z*z)
                        +2*gamma_1*gamma_2*(-1-6*z*z+15*z*z*z*z)
                        +gamma_1*gamma_1*(5-18*z*z+21*z*z*z*z)
                        +gamma_2*gamma_2*(5-18*z*z+21*z*z*z*z);
    double term_num = 24.*term_num_1*term_num_2;
    double term_denom = PI*PI*std::pow((-1*(gamma_1-gamma_2)*(gamma_1-gamma_2)*(-1+z*z)*(-1+z*z)
                        + std::pow((gamma_1+gamma_2+2*(lambda_e+lambda_Pi)-3*gamma_1*z*z-3*gamma_2*z*z),2)),7./2.);
    double res = -term_num/term_denom;
    return res/3.;
}

double analytic_dpi1dz_withoutm0(double z, double lambda_e, double lambda_Pi, double gamma_1, double gamma_2)
{
    double term_num_1 = -3.*std::pow((gamma_1-gamma_2),3.)*std::pow((-1.+z*z),4.);
    double term_num_2 = -12.*(gamma_1-gamma_2)*std::pow((-1+z*z),2)*
                          std::pow((gamma_1+gamma_2+2*(lambda_e+lambda_Pi) - 3*(gamma_1+gamma_2)*z*z),2);
    double term_num_3 = -2*(-1+3*z*z)*std::pow((gamma_1+gamma_2+2*(lambda_e+lambda_Pi) - 3*(gamma_1+gamma_2)*z*z),3);
    double term_num_4 = 3*std::pow((gamma_1-gamma_2),2)*std::pow((-1+z*z),2)*(-1+3*z*z)*(-2*(lambda_e+lambda_Pi)+(gamma_1+gamma_2)*(-1+3*z*z));
    double term_num = 4*(term_num_1 + term_num_2 + term_num_3 + term_num_4);
    double term_dom_1 = -std::pow((gamma_1-gamma_2),2)*std::pow((-1+z*z),2);
    double term_dom_2 = std::pow((gamma_1+gamma_2+2*(lambda_e+lambda_Pi)-3*(gamma_1+gamma_2)*z*z),2);
    double term_dom = PI*PI*std::pow(term_dom_1+term_dom_2,7./2.);
    double res = term_num/term_dom;
    return res;
}

double analytic_dpi2dz_withoutm0(double z, double lambda_e, double lambda_Pi, double gamma_1, double gamma_2)
{
    double term_num_1 = 3*std::pow((gamma_1-gamma_2),3)*std::pow((-1+z*z),4);
    double term_num_2 = 12*(gamma_1-gamma_2)*std::pow((-1+z*z),2)*
                          std::pow((gamma_1+gamma_2+2*(lambda_e+lambda_Pi) - 3*(gamma_1+gamma_2)*z*z),2);
    double term_num_3 = -2*(-1+3*z*z)*std::pow((gamma_1+gamma_2+2*(lambda_e+lambda_Pi) - 3*(gamma_1+gamma_2)*z*z),3);
    double term_num_4 = 3*std::pow((gamma_1-gamma_2),2)*std::pow((-1+z*z),2)*(-1+3*z*z)*(-2*(lambda_e+lambda_Pi)+(gamma_1+gamma_2)*(-1+3*z*z));
    double term_num = 4*(term_num_1 + term_num_2 + term_num_3 + term_num_4);
    double term_dom_1 = -std::pow((gamma_1-gamma_2),2)*std::pow((-1+z*z),2);
    double term_dom_2 = std::pow((gamma_1+gamma_2+2*(lambda_e+lambda_Pi)-3*(gamma_1+gamma_2)*z*z),2);
    double term_dom = PI*PI*std::pow(term_dom_1+term_dom_2,7./2.);
    double res = term_num/term_dom;
    return res;
}
