#ifndef __LJ
#define __LJ

#include "atoms.h"

double ljts(Atoms &atoms, double epsilon = 1.0, double sigma = 1.0,
            double cutoff = 5.0);

#endif // __LJ