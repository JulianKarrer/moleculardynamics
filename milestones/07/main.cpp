#include "atoms.h"
#include "verlet.h"
#include "xyz.h"
#include <neighbors.h>
#include <ducastelle.h>
#include <iostream>
#include <thermostat.h>

const size_t PRINT_EVERY_N {200};
const double TEMPERATURE{300};

void run_at_dt(double dt, double t_total, double thermos_until_t){
    // initialiye gold cluster from file
    auto [names, positions, velocities]{read_xyz_with_velocities("cluster_923.xyz")};
    Atoms atoms{Atoms(names, positions)};
    std::ofstream traj("traj.xyz");
    // create neighbour list and initialize forces
    NeighborList neighbor_list;
    neighbor_list.update(atoms, 10.0);
    (void) ducastelle(atoms, neighbor_list);
    double t{0};
    size_t i{0};
    while (t < t_total){
        verlet_step1(atoms.positions, atoms.velocities, atoms.forces, dt,
                     atoms.masses);
        neighbor_list.update(atoms, 10.0);
        double potential{ducastelle(atoms, neighbor_list)};
        verlet_step2(atoms.velocities, atoms.forces, dt, atoms.masses);
        if (t<thermos_until_t){
            berendsen_thermostat(atoms, TEMPERATURE, dt, thermos_until_t);
        }
        
        // compute total energy and print
        if (i%PRINT_EVERY_N==0){
            double kinetic = atoms.kinetic_energy();
            std::cout << "E_pot" << potential << "\tE_kin:" << kinetic
                << "\tSum:" << kinetic + potential << "\n";
            write_xyz(traj, atoms);
        }
        t += dt;
        i++;
    }
}

int main(int argc, char *argv[]) {
    double dt = 1.0;
    double t_total = dt*10000.;
    double thermos_until_t = dt*1000.;
    run_at_dt(dt, t_total, thermos_until_t);
   return 0;
}