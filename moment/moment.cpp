#include<iostream>
#include<cmath>
#include<vector>


double simpson_mom(double f_a, double f_b, double f_c, double h) {
    return (h / 6) * (f_a + 4 * f_c + f_b);
}

double adaptive_simpson_4phi(double (*func)(double, double, double, double, int, int,
                                        int, int, double, double), // 被积函数指针
                        double a, double b, double epsilon, int max_depth, int depth,
                        double lambda_e, double lambda_Pi, double gamma_1, double gamma_2,
                        int n, int n1, int n2, int n3, double z) {
    double c = (a + b) / 2.0;
    double h = b - a;

    // 计算被积函数的值
    double f_a = func(lambda_e, lambda_Pi, gamma_1, gamma_2,
                    n, n1, n2, n3, z, a);
    double f_b = func(lambda_e, lambda_Pi, gamma_1, gamma_2,
                    n, n1, n2, n3, z, b);
    double f_c = func(lambda_e, lambda_Pi, gamma_1, gamma_2,
                    n, n1, n2, n3, z, c);
    // 粗略积分
    double I1 = simpson_mom(f_a, f_b, f_c, h);

    // 分割区间后的积分
    double d = (a + c) / 2.0;
    double e = (c + b) / 2.0;
    double f_d = func(lambda_e, lambda_Pi, gamma_1, gamma_2,
                    n, n1, n2, n3, z, d);
    double f_e = func(lambda_e, lambda_Pi, gamma_1, gamma_2,
                    n, n1, n2, n3, z, e);
    double I2 = simpson_mom(f_a, f_c, f_d, h / 2) + simpson_mom(f_c, f_b, f_e, h / 2);
    if (std::abs(I2 - I1) < 15 * epsilon || depth >= max_depth) {
        return I2 + (I2 - I1) / 15.0; // 校正误差
    }
    //std::cout << depth << std::endl;
    //std::cout << a << std::endl;
    // 若精度不满足，递归细分
    return adaptive_simpson_4phi(func, a, c, epsilon/2, max_depth, depth+1,
                        lambda_e, lambda_Pi, gamma_1, gamma_2,n, n1, n2, n3, z) +
           adaptive_simpson_4phi(func, c, b, epsilon/2, max_depth, depth+1,
                        lambda_e, lambda_Pi, gamma_1, gamma_2,n, n1, n2, n3, z);
}



double adaptive_simpson_4z(double (*func)(double, double, double, double, int, int,
                                        int, int, double), // 被积函数指针
                        double a, double b, double epsilon, int max_depth, int depth,
                        double lambda_e, double lambda_Pi, double gamma_1, double gamma_2,
                        int n, int n1, int n2, int n3) {
    double c = (a + b) / 2.0;
    double h = b - a;

    // 计算被积函数的值
    double f_a = func(lambda_e, lambda_Pi, gamma_1, gamma_2,
                    n, n1, n2, n3, a);
    double f_b = func(lambda_e, lambda_Pi, gamma_1, gamma_2,
                    n, n1, n2, n3, b);
    double f_c = func(lambda_e, lambda_Pi, gamma_1, gamma_2,
                    n, n1, n2, n3, c);
    // 粗略积分
    double I1 = simpson_mom(f_a, f_b, f_c, h);

    // 分割区间后的积分
    double d = (a + c) / 2.0;
    double e = (c + b) / 2.0;
    double f_d = func(lambda_e, lambda_Pi, gamma_1, gamma_2,
                    n, n1, n2, n3, d);
    double f_e = func(lambda_e, lambda_Pi, gamma_1, gamma_2,
                    n, n1, n2, n3, e);
    double I2 = simpson_mom(f_a, f_c, f_d, h / 2) + simpson_mom(f_c, f_b, f_e, h / 2);
    if (std::abs(I2 - I1) < 15 * epsilon || depth >= max_depth) {
        return I2 + (I2 - I1) / 15.0; // 校正误差
    }
    //std::cout << depth << std::endl;
    //std::cout << a << std::endl;
    // 若精度不满足，递归细分
    return adaptive_simpson_4z(func, a, c, epsilon/2, max_depth, depth+1,
                        lambda_e, lambda_Pi, gamma_1, gamma_2,n, n1, n2, n3) +
           adaptive_simpson_4z(func, c, b, epsilon/2, max_depth, depth+1,
                        lambda_e, lambda_Pi, gamma_1, gamma_2,n, n1, n2, n3);
}




double dqqqdzdphi(double lambda_e, double lambda_pi, double gamma_1, double gamma_2,
                    int n, int n1, int n2, int n3, double z, double phi)
{
    double num_1 = std::tgamma(2+n+n1+n2+n3);
    double num_2 = std::pow(1.-z*z,(n1+n2)/2.)*std::pow(z,n3)*std::pow(std::cos(phi),n1)*std::pow(std::sin(phi),n2);
    double dom_ = lambda_e + lambda_pi - (gamma_1 + gamma_2)*z*z + (1.-z*z)*(gamma_1*std::pow(std::cos(phi),2) + gamma_2*std::pow(std::sin(phi),2));
    double dom = std::pow(dom_,2+n+n1+n2+n3);
    return num_1*num_2/dom;
}

double dqqqdz(double lambda_e, double lambda_pi, double gamma_1, double gamma_2,
                    int n, int n1, int n2, int n3, double z)

{    
    const double a = 0.0; 
    const double b = 2.0 * 3.1415926; 
    double res = adaptive_simpson_4phi(dqqqdzdphi, a, b, 1e-7, 30000, 0,
                        lambda_e, lambda_pi, gamma_1, gamma_2,
                        n, n1, n2, n3, z );
    return res;
}


double qqq(double lambda_e, double lambda_pi, double gamma_1, double gamma_2,
              int n, int n1, int n2, int n3) 
{
    const double a = -1.0; // 积分下限
    const double b = 1.0; // 积分上限
    double res = adaptive_simpson_4z(dqqqdz, a, b, 1e-7, 30000, 0,
                        lambda_e, lambda_pi, gamma_1, gamma_2,
                        n, n1, n2, n3);
    double prefactor = std::pow((2*3.1415926*0.1973),3);
    return res/prefactor;
                        
}

double qqq0(double T, int n, int n1, int n2, int n3)
{
    double num_1 = std::tgamma(2+n+n1+n2+n3)*(1+std::pow(-1,n1))*(1+std::pow(-1,n2))*(1+std::pow(-1,n3));
    double num_2 = std::tgamma((1.+n1)/2.)*std::tgamma((1.+n2)/2.)*std::tgamma((1.+n3)/2.);
    double dom = 4.*std::tgamma((3.+n1+n2+n3)/2.)*std::pow((1./T),(2+n+n1+n2+n3));
    double prefactor = std::pow((2*3.1415926*0.1973),3);
    return num_1*num_2/(dom*prefactor);
}


int main(int argc, char* argv[]) {

    const double lambda_e = atof(argv[1]);
    const double lambda_pi= atof(argv[2]);
    const double gamma_1= atof(argv[3]);
    const double gamma_2= atof(argv[4]);
    const double beta = atof(argv[5]);
    const int n = atof(argv[6]);
    const int n1 = atof(argv[7]);
    const int n2 = atof(argv[8]);
    const int n3 = atof(argv[9]);


    double moment_fh = qqq(lambda_e, lambda_pi, gamma_1, gamma_2, n, n1, n2, n3);
    std::cout << "moment_fh: " << moment_fh << std::endl;
    double aa = qqq0(1./beta, n, n1, n2, n3);
    std::cout << "moment_f0: " << aa << std::endl;
    return 0;
}