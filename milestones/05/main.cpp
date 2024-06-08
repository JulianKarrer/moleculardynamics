#include "atoms.h"
#include "lj_direct_summation.h"
#include "verlet.h"
#include "xyz.h"
#include <cmath>
#include <iostream>
#include <thermostat.h>

const double DT{0.001};
const double SIGMA{1.44};
const double TEMPERATURE{300.};

int main(int argc, char *argv[]) {
    // initialize lattice and trajectory file
    Atoms atoms{Atoms(27, SIGMA * pow(2, 1. / 6.))};
    std::ofstream traj("lattice.xyz");
    for (int i{0}; i < 10000; ++i) {
        double potential_energy = lj_direct_summation(atoms, 1.0, SIGMA);
        verlet_step1(atoms.positions, atoms.velocities, atoms.forces, DT,
                     atoms.masses);
        potential_energy += lj_direct_summation(atoms, 1.0, SIGMA);
        verlet_step2(atoms.velocities, atoms.forces, DT, atoms.masses);
        berendsen_thermostat(atoms, TEMPERATURE, DT, 100 * DT);

        std::cout << potential_energy << " kin:" << atoms.kinetic_energy()
                  << "\n";
        if (i % 100 == 0) {
            // write every 1000th simulation step to disk
            write_xyz(traj, atoms);
        }
    }
    write_xyz(traj, atoms);
    traj.close();
    return 0;
}
