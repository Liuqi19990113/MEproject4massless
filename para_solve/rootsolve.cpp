#include<iostream>
#include<random>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>
#include"elements.h"
#include"Integral.h"
#include"rootsolve.h"



int func(const gsl_vector * x, void *params,
	gsl_vector * f)
{
    double prefactor = 2;
	const double e_ = ((struct rparams *) params)->e_real; 
 	const double pi1_ = ((struct rparams *) params)->pi1_real; 
	const double pi2_ = ((struct rparams *) params)->pi2_real;
     

	const double lambda_e = gsl_vector_get(x, 0); 
    const double gamma_1 =  gsl_vector_get(x, 1);
    const double gamma_2 = gsl_vector_get(x, 2);

    const double e = (prefactor*integral_z_for_e(lambda_e, 0., gamma_1, gamma_2) - e_)/e_;
    const double pi1 = (prefactor*integral_z_for_pi1(lambda_e, 0., gamma_1, gamma_2) - pi1_)/pi1_;
    const double pi2 = (prefactor*integral_z_for_pi2(lambda_e, 0., gamma_1, gamma_2) - pi2_)/pi2_;
    //std::cout << prefactor*integral_z_for_e(lambda_e, 0., gamma_1, gamma_2)
    //            << " " << integral_z_for_pi1(lambda_e, 0., gamma_1, gamma_2) << " "
    //            << integral_z_for_pi1(lambda_e, 0., gamma_1, gamma_2) << std::endl;
    //std::cout << e << " " << pi1 << " " << pi2 << std::endl;

	gsl_vector_set(f, 0, e);
	gsl_vector_set(f, 1, pi1);
	gsl_vector_set(f, 2, pi2);
    return GSL_SUCCESS;
}

int print_state_f(size_t iter, gsl_multiroot_fsolver * s)
{
	//printf("iter = %3u Lambda, gamma1, gamma2 = % .8f  % .8f % .8f"
	//	"  f(x) = % .5e % .5e % .5e \n",
	//	iter,
      printf("% .8f  % .8f % .8f" " % .5e % .5e % .5e \n",
		gsl_vector_get(s->x, 0),
		gsl_vector_get(s->x, 1),
		gsl_vector_get(s->x, 2),
		gsl_vector_get(s->f, 0),
		gsl_vector_get(s->f, 1),
        gsl_vector_get(s->f, 2));
	return 0;
}

int solve(struct rparams p, double x_init[])
{
    const gsl_multiroot_fsolver_type *T;
    gsl_multiroot_fsolver *s;
	int status;
	size_t i, iter = 0;

	const size_t n = 3;

    gsl_multiroot_function f = {&func,n, &p };

    gsl_vector *x = gsl_vector_alloc (n);

    gsl_vector_set (x, 0, x_init[0]);
    gsl_vector_set (x, 1, x_init[1]);
    gsl_vector_set (x, 2, x_init[2]);


    T = gsl_multiroot_fsolver_hybrids;
    s = gsl_multiroot_fsolver_alloc (T, n);
    gsl_multiroot_fsolver_set (s, &f, x);

    //print_state_f (iter, s);

    do
    {
      iter++;
      status = gsl_multiroot_fsolver_iterate (s);

      //print_state_f (iter, s);

      if (status)   /* check if solver is stuck */
        break;

        status =
        gsl_multiroot_test_residual (s->f, 1e-6);
    }while (status == GSL_CONTINUE && iter < 1000);

    //printf ("status = %s\n", gsl_strerror (status));
    //printf (gsl_strerror (status));
    print_state_f (iter, s);
    // std::cout << std::fixed << std::setprecision(7) << " "
    //           << gsl_vector_get(s->x, 0) << "  "
	// 	      <<gsl_vector_get(s->x, 1)  << "  "
	// 	      <<gsl_vector_get(s->x, 2) << "  "
	// 	      <<gsl_vector_get(s->x, 3);


    gsl_multiroot_fsolver_free (s);
    gsl_vector_free (x);
    return 0;
}