#include "atoms.h"
#include "lj.h"
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

int main(int argc, char *argv[]) {
    // initialize lattice and trajectory file
    Atoms atoms{Atoms(1000, SIGMA * pow(2, 1. / 6.))};
    std::ofstream traj("lattice.xyz");

    (void)ljts(atoms, 1.0, SIGMA, CUTOFF);
    for (int i{0}; i < 5000; ++i) {
        // compute lj forces and integrate with velocity verlet
        verlet_step1(atoms.positions, atoms.velocities, atoms.forces, DT,
                     atoms.masses);
        double potential_energy{ljts(atoms, 1.0, SIGMA, CUTOFF)};
        verlet_step2(atoms.velocities, atoms.forces, DT, atoms.masses);
        // quilibrate the system with a berendesen thermostat
        berendsen_thermostat(atoms, TEMPERATURE, DT, 100 * DT);

        std::cout << "E_pot:" << potential_energy
                  << "\tE_kin:" << atoms.kinetic_energy() << "\n";
        if (i % 100 == 0) {
            // write every 100th simulation step to disk
            write_xyz(traj, atoms);
        }
    }
    write_xyz(traj, atoms);
    traj.close();
    return 0;
}
