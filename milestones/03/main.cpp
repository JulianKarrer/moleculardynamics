#include "atoms.h"
#include "lj_direct_summation.h"
#include "verlet.h"
#include "xyz.h"
#include <cmath>
#include <iostream>

inline void print_energy(double potential, double kinetic) {
    std::cout << potential + kinetic << "\n";
}

int main(int argc, char *argv[]) {
    auto [names, positions, velocities]{read_xyz_with_velocities("lj54.xyz")};
    Atoms atoms{Atoms(positions, velocities)};
    double dt{0.001};
    double m{1.};
    std::ofstream traj("traj.xyz");
    for (int i{0}; i < 100000; ++i) {
        verlet_step1(atoms.positions, atoms.velocities, atoms.forces, dt, m);
        double potential_energy = lj_direct_summation(atoms);
        verlet_step2(atoms.velocities, atoms.forces, dt, m);
        potential_energy += lj_direct_summation(atoms);
        print_energy(potential_energy, atoms.kinetic_energy(m));
        if (i % 1000 == 0) {
            // write every 1000th simulation step to disk
            write_xyz(traj, atoms);
        }
    }
    traj.close();
    return 0;
}
