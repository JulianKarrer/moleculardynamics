#include "atoms.h"
#include "ljts.h"
#include "verlet.h"
#include "xyz.h"
#include <cmath>
#include <iostream>
#include <neighbors.h>
#include <thermostat.h>

const double DT{0.001};
const double SIGMA{0.8};
const double CUTOFF{5. * SIGMA};
const double TEMPERATURE{300.};

/// run simulation for given relaxation time and push temperature at each
/// timestep into the given temperatures vector
void temperature_graph(std::vector<double> &temps, double relaxation) {
    Atoms atoms{Atoms(100, SIGMA * pow(2, 1. / 6.))};
    (void)ljts(atoms, 1.0, SIGMA, CUTOFF);
    for (int i{0}; i < 1000; i++) {
        verlet_step1(atoms.positions, atoms.velocities, atoms.forces, DT,
                     atoms.masses);
        (void)ljts(atoms, 1.0, SIGMA, CUTOFF);
        verlet_step2(atoms.velocities, atoms.forces, DT, atoms.masses);
        berendsen_thermostat(atoms, TEMPERATURE, DT, relaxation);
        temps.push_back(temperature_cur(atoms));
    }
}

int main(int argc, char *argv[]) {
    // EQUILIBRATE A LATTICE AND SEE THE RESULTING TRAJECTORY ~~~~~~~~~~~
    // initialize lattice and trajectory file
    Atoms atoms{Atoms(1000, SIGMA * pow(2, 1. / 6.))};
    std::ofstream traj("lattice.xyz");

    (void)ljts(atoms, 1.0, SIGMA, CUTOFF);
    for (int i{0}; i < 5000; i++) {
        // compute lj forces and integrate with velocity verlet
        verlet_step1(atoms.positions, atoms.velocities, atoms.forces, DT,
                     atoms.masses);
        double potential_energy{ljts(atoms, 1.0, SIGMA, CUTOFF)};
        verlet_step2(atoms.velocities, atoms.forces, DT, atoms.masses);
        // eyquilibrate the system with a berendesen thermostat
        berendsen_thermostat(atoms, TEMPERATURE, DT, 100 * DT);

        // write every 100th simulation step to disk
        if (i % 100 == 0) {
            std::cout << "E_pot:" << potential_energy
                      << "\tE_kin:" << atoms.kinetic_energy() << "\n";
            write_xyz(traj, atoms);
        }
    }
    write_xyz(traj, atoms);
    traj.close();

    // PLOT RELAXATION TIMES AND TEMPERATURES ~~~~~~~~~~~~~~~~~~~~~~~~~~~

    std::ofstream file("thermostat_temps.csv");
    file << "relaxation time,current_time,temperature" << std::endl;

    for (size_t factor{10}; factor < 100; factor += 10) {
        std::cout << factor << std::endl;
        // run simulation with given relaxation time
        double relaxation{(double)factor * DT};
        std::vector<double> temps;
        temperature_graph(temps, relaxation);
        // dump the resulting temperatures per timestep into a csv
        for (size_t current_time{0}; current_time < temps.size();
             current_time++) {
            file << relaxation << "," << (double)(current_time + 1) * DT << ","
                 << temps[current_time] << std::endl;
        }
    }

    return 0;
}
