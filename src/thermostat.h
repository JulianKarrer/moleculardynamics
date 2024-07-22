#ifndef __THERMOSTAT
#define __THERMOSTAT

#include "atoms.h"

/// The exact Boltzman constant in units of eV/K to 9 decimal places.
const double KB_EV_K{8.617333262e-5};
/// 3/2 times the exact Boltzman constant in units of eV/K to 9 decimal places.
const double KB_EV_K_3_2{1.5 * KB_EV_K};

/// @brief A berendsen thermostat. Directly operates on the velocities of the
/// `Atoms` passed in.
/// @param atoms atoms to modify
/// @param temperature target temperature
/// @param dt timestep
/// @param relaxation_time relaxation time, on the order of tau >> dt
void berendsen_thermostat(Atoms &atoms, double temperature, double dt,
                          double relaxation_time);

/// @brief Like the berendsen thermostat, but used on a decomposed set of
/// `Atoms` as used in the MPI parallelization
/// @param atoms atoms to modify
/// @param temperature target temperature
/// @param dt timestep
/// @param relaxation_time relaxation time, on the order of tau >> dt
/// @param nb_local first `nb_local` atoms in `Atoms` are assumed not to be
/// ghost atoms
/// @param mass the mass of each atom, which assumes they have the same mass
void berendsen_thermostat_decomposed(Atoms &atoms, double temperature,
                                     double dt, double relaxation_time,
                                     size_t nb_local, double mass);

/// @brief Current temperature of the system.
/// @param atoms `Atoms` to measure the temperature of
/// @return current temperature in Kelvins
double temperature_cur(const Atoms atoms);

#endif // __THERMOSTAT