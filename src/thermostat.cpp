#include "thermostat.h"

double temperature_cur(Atoms &atoms) {
    return atoms.kinetic_energy() / (KB_EV_K_3_2 * atoms.nb_atoms());
}

void berendsen_thermostat(Atoms &atoms, double temperature, double dt,
                          double relaxation_time) {
    double lambda{sqrt(1 + (temperature / temperature_cur(atoms) - 1) * dt /
                               relaxation_time)};
    atoms.velocities *= lambda;
}