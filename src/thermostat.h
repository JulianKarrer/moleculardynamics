#ifndef __THERMOSTAT
#define __THERMOSTAT

#include "atoms.h"

/// The exact Boltzman constant in units of eV/K to 9 decimal places.
const double KB_EV_K{8.617333262e-5};
/// 3/2 times the exact Boltzman constant in units of eV/K to 9 decimal places.
const double KB_EV_K_3_2{1.5 * KB_EV_K};

void berendsen_thermostat(Atoms &atoms, double temperature, double dt,
                          double relaxation_time);

double temperature_cur(Atoms &atoms);

#endif // __THERMOSTAT