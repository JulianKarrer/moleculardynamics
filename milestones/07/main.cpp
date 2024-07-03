#include "atoms.h"
#include "verlet.h"
#include "xyz.h"
#include <ducastelle.h>
#include <iomanip>
#include <iostream>
#include <neighbors.h>
#include <thermostat.h>

/// @brief Perform one simulation step
/// @param atoms the atoms of the system to propagate
/// @param dt timestep size in fs
/// @param neighbour_list some neighbour list that will be updated
/// @returns the current potential energy
double step(Atoms &atoms, double dt, NeighborList &neighbour_list) {
    verlet_step1(atoms.positions, atoms.velocities, atoms.forces, dt,
                 atoms.masses);
    neighbour_list.update(atoms, 10.0);
    double potential{ducastelle(atoms, neighbour_list)};
    verlet_step2(atoms.velocities, atoms.forces, dt, atoms.masses);
    return potential;
};

/// @brief Run a simulation for a given timestep, tracking the hamiltonian and
/// writing it to a csv file for analysis
/// @param dt the timestep to use
/// @param total_time the total time to simulate the system
/// @param dt_eam_file the file to write the results to
void run_at_dt(double dt, double total_time, std::ofstream &dt_eam_file) {
    // initialiye gold cluster from file
    auto [names, positions,
          velocities]{read_xyz_with_velocities("cluster_923.xyz")};
    Atoms atoms{Atoms(names, positions)};
    // create neighbour list and initialize forces
    NeighborList neighbor_list;
    neighbor_list.update(atoms, 10.0);
    (void)ducastelle(atoms, neighbor_list);
    double t{0};
    while (t < total_time) {
        double potential{step(atoms, dt, neighbor_list)};
        double kinetic = atoms.kinetic_energy();
        dt_eam_file << std::fixed << std::setprecision(2) << dt << ","
                    << std::fixed << std::setprecision(10) << t << ","
                    << kinetic + potential << "," << potential << "," << kinetic
                    << std::endl;
        t += dt;
    }
}

void debug_print(Atoms atoms, double e_pot) {
    double e_kin{atoms.kinetic_energy()};
    std::cout << std::fixed << std::setprecision(5)
              << "System heated up to: " << temperature_cur(atoms)
              << "K | E_total: " << e_pot + e_kin << "eV | E_pot: " << e_pot
              << "eV | E_kin: " << e_kin << "eV" << std::endl;
}

/// @brief Specify one or more arguments to run a test on what timestep size
/// is appropriate for this potential. Otherwise, run the simulation and
/// collect data.
int main(int argc, char *argv[]) {
    if (argc > 1) {
        // FIND OPTIMAL TIMESTEP SIZE

        // total time to simulate the system for
        double total_time{1000.};
        // open a csv file for the results, write the header
        std::ofstream dt_eam_file("dt_eam.csv");
        dt_eam_file << "dt,t,Hamiltonian,E_pot,E_kin," << std::endl;
        // try some timestep sizes, record the time evolution of the hamiltonian
        // for each
        for (double dt : {0.25, 0.5, 1.0, 5.0, 10.}) {
            run_at_dt(dt, total_time, dt_eam_file);
        }
    }

    // EAM ISOCAHEDRON ENERGY-TEMPERATURE GRAPH
    // settings
    double dt{0.1};
    size_t tau_relax{5000}; // in units of dt
    double from_temp{500};
    double until_temp{1200};
    double increase_q_by{0.5};
    std::ofstream gold_t_e("gold_t_e.csv");

    // initialize gold cluster from file
    auto [names, positions,
          velocities]{read_xyz_with_velocities("cluster_923.xyz")};
    Atoms atoms{Atoms(names, positions)};

    // create neighbour list and initialize forces
    NeighborList neighbour_list;
    neighbour_list.update(atoms, 10.0);
    (void)ducastelle(atoms, neighbour_list);
    double current_temp{0.};

    // first, heat up the system gently to some initial temperature so we don't
    // have to start from zero
    double e_pot{0.0};
    for (size_t i{0}; i < 10 * tau_relax; i++) {
        // perform timestep
        e_pot = step(atoms, dt, neighbour_list);
        berendsen_thermostat(atoms, from_temp, dt, tau_relax * dt);
        if (i % 500 == 0) {
            debug_print(atoms, e_pot);
        }
    }
    debug_print(atoms, e_pot);
    std::cout << "MEASUREMENTS STARTING" << std::endl;

    // then, repeatedly increase temperature, wait for tau_relax, measure energy
    // and temperature for tau_relax
    gold_t_e << "n,e_total,temperature,e_pot,e_kin" << std::endl;
    // might as well record the trajectory to check plausibility of results
    std::ofstream traj("gold_t_e_traj.xyz");
    while (current_temp < until_temp) {
        // bump up energy
        atoms.increase_kinetic_energy(increase_q_by);
        // wait for system to relax
        for (size_t i{0}; i <= tau_relax; i++) {
            (void)step(atoms, dt, neighbour_list);
        }
        // measure T and E
        double avg_e_tot{0.0};
        double avg_e_pot{0.0};
        double avg_e_kin{0.0};
        double avg_t{0.0};
        for (size_t i{0}; i <= tau_relax; i++) {
            e_pot = step(atoms, dt, neighbour_list);
            double e_kin{atoms.kinetic_energy()};
            avg_e_tot += e_pot + e_kin;
            avg_e_pot += e_pot;
            avg_e_kin += e_kin;
            avg_t += temperature_cur(atoms);
        }
        // write data to csv
        double normalize{1. / ((double)tau_relax)};
        gold_t_e << atoms.nb_atoms() << "," << std::fixed
                 << std::setprecision(15) << avg_e_tot * normalize << ","
                 << avg_t * normalize << "," << avg_e_pot * normalize << ","
                 << avg_e_kin * normalize << std::endl;
        write_xyz(traj, atoms);
        // update current temperature
        current_temp = temperature_cur(atoms);
        // print to std:out to watch the process
        std::cout << std::fixed << std::setprecision(5)
                  << "System heated up to: " << avg_t * normalize
                  << "K | E_total: " << avg_e_tot * normalize
                  << "eV | E_pot: " << avg_e_pot * normalize
                  << "eV | E_kin: " << avg_e_kin * normalize << "eV"
                  << std::endl;
    }

    return 0;
}