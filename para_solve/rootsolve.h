#ifndef ROOTSOLVE_H
#define ROOTSOLVE_H
#include<random>
#include<iostream>
#include<iomanip>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>
#include"elements.h"
#include"Integral.h"

struct rparams
{
	double e_real;
    double pi1_real;
    double pi2_real;
};

int func(const gsl_vector * x, void *params,
	gsl_vector * f);

int print_state_f(size_t iter, gsl_multiroot_fsolver * s);

int solve(struct rparams p, double x_init[]);

#endif