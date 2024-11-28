#include<iostream>
#include<cmath>
#include<iomanip>
#include"elements.h"
#include"Integral.h"
#include"rootsolve.h"

int main(int argc, char* argv[])
{
    //  double e =integral_z_for_e(6.9564, 0., -1.0927, 2.1855);
    //  double p1 = integral_z_for_pi1(6.9564, 0., -1.0927, 2.1855);
    //  double p2 = integral_z_for_pi2(6.9564, 0., -1.0927, 2.1855);
    //  std::cout << e << " " << p1 << " " << p2 << std::endl;


    //  double e_ =integral_z_for_e(9.1582, 0., -1.4461, 2.8807);
    //  double p1_ = integral_z_for_pi1(9.1582, 0., -1.4461, 2.8807);
    //  double p2_ = integral_z_for_pi2(9.1582, 0., -1.4461, 2.8807);
    //  std::cout << e_ << " " << p1_ << " " << p2_ << std::endl;



    //  std::cout << (e-atof(argv[1]))/atof(argv[1]) << " " 
    //  << (p1 - atof(argv[2]))/atof(argv[2]) << " " 
    //  << (p2 - atof(argv[3]))/atof(argv[3]) << std::endl;
    const double e_real = atof(argv[1]);
    const double pi1_real = atof(argv[2]);
    const double pi2_real = atof(argv[3]);
    double x_guess[4] = {atof(argv[4]),atof(argv[5]),atof(argv[6])};
    struct rparams real_value = {e_real, pi1_real, pi2_real};
    solve(real_value, x_guess);

}