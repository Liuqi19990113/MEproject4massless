#include<cmath>
#include"Integral.h"

double integral_q(double (*f)(double, double, double, double, double, double, double), double m0, double q_max, double z, double lambda_e, double lambda_Pi, double gamma_1, double gamma_2)
{
    double dq = 0.2;
    double integral_q_res = 0;
    long long n = q_max/dq;
    double fa = f(m0 ,z, lambda_e, lambda_Pi, gamma_1, gamma_2, 0);
    double fb = f(m0, z, lambda_e, lambda_Pi, gamma_1, gamma_2, q_max);
    integral_q_res += fa;
    integral_q_res += fb;

    //simpson's rule
    long long j = 0;
    for(long long i = 1; i < n; i++)
    {
        double q_i = i*dq;
        double q_imhalf = dq*double(i+j)/2;
        j = i;
        double f_i = f(m0, z, lambda_e, lambda_Pi, gamma_1, gamma_2, q_i);
        double f_imhalf = f(m0, z, lambda_e, lambda_Pi, gamma_1, gamma_2, q_imhalf);
        integral_q_res += (4*f_imhalf);
        integral_q_res += (2*f_i);
    }

    integral_q_res += 4*f(m0, z, lambda_e, lambda_Pi, gamma_1, gamma_2, q_max-dq/2);
    double h = q_max/n;
    return integral_q_res*h/6;
}


double simpson(double f_a, double f_b, double f_c, double h) {
    return (h / 6) * (f_a + 4 * f_c + f_b);
}

double adaptive_simpson(double (*func)(double, double, double, double, double), // 被积函数指针
                        double a, double b, double epsilon, int max_depth, int depth,
                        double lambda_e, double lambda_Pi, double gamma_1, double gamma_2) {
    double c = (a + b) / 2.0;
    double h = b - a;

    // 计算被积函数的值
    double f_a = func(a, lambda_e, lambda_Pi, gamma_1, gamma_2);
    double f_b = func(b, lambda_e, lambda_Pi, gamma_1, gamma_2);
    double f_c = func(c, lambda_e, lambda_Pi, gamma_1, gamma_2);

    // 粗略积分
    double I1 = simpson(f_a, f_b, f_c, h);

    // 分割区间后的积分
    double d = (a + c) / 2.0;
    double e = (c + b) / 2.0;
    double f_d = func(d, lambda_e, lambda_Pi, gamma_1, gamma_2);
    double f_e = func(e, lambda_e, lambda_Pi, gamma_1, gamma_2);
    double I2 = simpson(f_a, f_c, f_d, h / 2) + simpson(f_c, f_b, f_e, h / 2);
    if (std::abs(I2 - I1) < 15 * epsilon || depth >= max_depth) {
        return I2 + (I2 - I1) / 15.0; // 校正误差
    }
    //std::cout << depth << std::endl;
    //std::cout << a << std::endl;
    // 若精度不满足，递归细分
    return adaptive_simpson(func, a, c, epsilon / 2, max_depth, depth + 1, lambda_e, lambda_Pi, gamma_1, gamma_2) +
           adaptive_simpson(func, c, b, epsilon / 2, max_depth, depth + 1, lambda_e, lambda_Pi, gamma_1, gamma_2);
}


double integral_z_for_e(double lambda_e, double lambda_Pi, double gamma_1, double gamma_2)
{
    // double z_max = 1;
    // double dz = 0.01;
    // long long m = z_max/dz;
    // //double m0 = 1.;
    // //double q_max = 1000.*m0;
    // double integral_z_res = 0;
    // double fa = analytic_dedz_withoutm0(0., lambda_e, lambda_Pi, gamma_1, gamma_2);
    // //std::cout << integral_q(dedqdz_withm0, m0, q_max , 0., lambda_e, lambda_Pi, gamma_1, gamma_2) << std::endl;
    // //std::cout << integral_q(dedqdz_withm0, 0., q_max , 0., lambda_e, lambda_Pi, gamma_1, gamma_2) << std::endl;
    // //std::cout << analytic_dedz_withoutm0(0., lambda_e, lambda_Pi, gamma_1, gamma_2) << std::endl;
    // //std::cout << fa <<std::endl;
    // double fb = analytic_dedz_withoutm0(z_max, lambda_e, lambda_Pi, gamma_1, gamma_2);
    // //std::cout << fb <<std::endl;    
    // integral_z_res += fa;
    // integral_z_res += fb;
    // //simpson's rule
    // long long j = 0;
    // for(long long i = 1; i < m; i++)
    // {
    //     double z_i = i*dz;
    //     double z_imhalf = dz*double(i+j)/2;
    //     j = i;
    //     double f_i = analytic_dedz_withoutm0(z_i, lambda_e, lambda_Pi, gamma_1, gamma_2);
    //     double f_imhalf =  analytic_dedz_withoutm0(z_imhalf, lambda_e, lambda_Pi, gamma_1, gamma_2);
    //     //std::cout << f_imhalf << std::endl;
    //     //if (i == 20){std::cout << f_i << std::endl;}
    //     integral_z_res += (4*f_imhalf);
    //     integral_z_res += (2*f_i);
    //     //std::cout << integral_z_res << std::endl;
    // }
    // integral_z_res += 4*(analytic_dedz_withoutm0(z_max-dz/2, lambda_e, lambda_Pi, gamma_1, gamma_2));
    // double h = z_max/m;
    // //{std::cout << integral_z_res << std::endl;}
    // return integral_z_res*h/(6*(std::pow(0.197,3)));
    double result = adaptive_simpson(analytic_dedz_withoutm0, 0., 1., 1e-8, 100000, 0, lambda_e, lambda_Pi, gamma_1, gamma_2);
    return result/(std::pow(0.1973,3));
}

double integral_z_for_P(double lambda_e, double lambda_Pi, double gamma_1, double gamma_2)
{
    // double z_max = 1;
    // double dz = 0.01;
    // long long m = z_max/dz;
    // //double m0 = 1.;
    // //double q_max = 1000.*m0;
    // double integral_z_res = 0;
    // double fa = analytic_dPdz_withoutm0(0., lambda_e, lambda_Pi, gamma_1, gamma_2);
    // //std::cout << integral_q(dPdqdz_withm0, m0, q_max , 0., lambda_e, lambda_Pi, gamma_1, gamma_2) << std::endl;
    // //std::cout << integral_q(dPdqdz_withm0, 0., q_max , 0., lambda_e, lambda_Pi, gamma_1, gamma_2) << std::endl;
    // //std::cout << analytic_dPdz_withoutm0(0., lambda_e, lambda_Pi, gamma_1, gamma_2) << std::endl;
    // //std::cout << fa <<std::endl;
    // double fb = analytic_dPdz_withoutm0(z_max, lambda_e, lambda_Pi, gamma_1, gamma_2);
    // integral_z_res += fa;
    // integral_z_res += fb;
    // //simpson's rule
    // long long j = 0;
    // for(long long i = 1; i < m; i++)
    // {
    //     double z_i = i*dz;
    //     double z_imhalf = dz*double(i+j)/2;
    //     j = i;
    //     double f_i =  analytic_dPdz_withoutm0(z_i, lambda_e, lambda_Pi, gamma_1, gamma_2);
    //     double f_imhalf = analytic_dPdz_withoutm0(z_imhalf, lambda_e, lambda_Pi, gamma_1, gamma_2);
    //     integral_z_res += (4*f_imhalf);
    //     integral_z_res += (2*f_i);
    // }
    // integral_z_res += 4*(analytic_dPdz_withoutm0(z_max-dz/2, lambda_e, lambda_Pi, gamma_1, gamma_2));
    // double h = z_max/m;
    // return integral_z_res*h/(6*(std::pow(0.197,3)));
    double result = adaptive_simpson(analytic_dPdz_withoutm0, 0, 1, 1e-8, 100000, 0, lambda_e, lambda_Pi, gamma_1, gamma_2);
    return result/(std::pow(0.1973,3));
}

double integral_z_for_pi1(double lambda_e, double lambda_Pi, double gamma_1, double gamma_2)
{
    // double z_max = 1;
    // double dz = 0.01;
    // long long m = z_max/dz;
    // //double m0 = 1.;
    // //double q_max = 1000.*m0;
    // double integral_z_res = 0;
    // double fa = analytic_dpi1dz_withoutm0(0., lambda_e, lambda_Pi, gamma_1, gamma_2);
    // //std::cout << integral_q(dpi1dqdz_withm0, 0., q_max , 0., lambda_e, lambda_Pi, gamma_1, gamma_2) << std::endl;
    // //std::cout << analytic_dpi1dz_withoutm0(0., lambda_e, lambda_Pi, gamma_1, gamma_2) << std::endl;
    // double fb = analytic_dpi1dz_withoutm0(z_max, lambda_e, lambda_Pi, gamma_1, gamma_2);
    // integral_z_res += fa;
    // integral_z_res += fb;
    // //simpson's rule
    // long long j = 0;
    // for(long long i = 1; i < m; i++)
    // {
    //     double z_i = i*dz;
    //     double z_imhalf = dz*double(i+j)/2;
    //     j = i;
    //     double f_i =  analytic_dpi1dz_withoutm0(z_i, lambda_e, lambda_Pi, gamma_1, gamma_2);
    //     double f_imhalf = analytic_dpi1dz_withoutm0(z_imhalf, lambda_e, lambda_Pi, gamma_1, gamma_2);
    //     integral_z_res += (4*f_imhalf);
    //     integral_z_res += (2*f_i);
    // }
    // integral_z_res += 4*(analytic_dpi1dz_withoutm0(z_max-dz/2, lambda_e, lambda_Pi, gamma_1, gamma_2));
    // double h = z_max/m;
    // return integral_z_res*h/(6*(std::pow(0.197,3)));
    double result = adaptive_simpson(analytic_dpi1dz_withoutm0, 0, 1, 1e-8, 100000, 0, lambda_e, lambda_Pi, gamma_1, gamma_2);
    return result/(std::pow(0.1973,3));
}

double integral_z_for_pi2(double lambda_e, double lambda_Pi, double gamma_1, double gamma_2)
{
    // double z_max = 1;
    // double dz = 0.01;
    // long long m = z_max/dz;
    // //double m0 = 1.;
    // //double q_max = 1000.*m0;
    // double integral_z_res = 0;
    // double fa = analytic_dpi2dz_withoutm0(0., lambda_e, lambda_Pi, gamma_1, gamma_2);
    // //std::cout << integral_q(dpi2dqdz_withm0, 0., q_max , 0., lambda_e, lambda_Pi, gamma_1, gamma_2) << std::endl;
    // //std::cout << analytic_dpi2dz_withoutm0(0., lambda_e, lambda_Pi, gamma_1, gamma_2) << std::endl;
    // double fb = analytic_dpi2dz_withoutm0(z_max, lambda_e, lambda_Pi, gamma_1, gamma_2);
    // integral_z_res += fa;
    // integral_z_res += fb;
    // //simpson's rule
    // long long j = 0;
    // for(long long i = 1; i < m; i++)
    // {
    //     double z_i = i*dz;
    //     double z_imhalf = dz*double(i+j)/2;
    //     j = i;
    //     double f_i =  analytic_dpi2dz_withoutm0(z_i, lambda_e, lambda_Pi, gamma_1, gamma_2);
    //     double f_imhalf = analytic_dpi2dz_withoutm0(z_imhalf, lambda_e, lambda_Pi, gamma_1, gamma_2);
    //     integral_z_res += (4*f_imhalf);
    //     integral_z_res += (2*f_i);
    // }
    // integral_z_res += 4*(analytic_dpi2dz_withoutm0(z_max-dz/2, lambda_e, lambda_Pi, gamma_1, gamma_2));
    // double h = z_max/m;
    // return integral_z_res*h/(6*(std::pow(0.197,3)));
    double result = adaptive_simpson(analytic_dpi2dz_withoutm0, 0, 1, 1e-8, 100000, 0, lambda_e, lambda_Pi, gamma_1, gamma_2);
    return result/(std::pow(0.1973,3));
}
