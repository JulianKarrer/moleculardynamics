#include "thermostat.h"

double temperature_cur(Atoms &atoms) {
    return atoms.kinetic_energy() / (KB_EV_K_3_2 * atoms.nb_atoms());
}

void berendsen_thermostat(Atoms &atoms, double temperature, double dt,
                          double relaxation_time) {
    double temp_cur{temperature_cur(atoms)};
    // avoid division by zero, implement berendsen thermostat as described in
    // milestone 5
    double lambda{temp_cur > 0. ? sqrt(1 + (temperature / temp_cur - 1) * dt /
                                               relaxation_time)
                                : 1.};
    atoms.velocities *= lambda;
}
