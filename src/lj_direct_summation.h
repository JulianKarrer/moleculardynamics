#ifndef __LJ_DIRECT_SUMMATION
#define __LJ_DIRECT_SUMMATION

#include "atoms.h"
double lj_potential(double r, double epsilon, double sigma);
double lj_potential_derivative(double r, double epsilon, double sigma);
double lj_direct_summation(Atoms &atoms, double epsilon = 1.0,
                           double sigma = 1.0);

#endif // __LJ_DIRECT_SUMMATION