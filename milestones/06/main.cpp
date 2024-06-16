#include "atoms.h"
#include "lj_direct_summation.h"
#include "ljts.h"
#include "verlet.h"
#include "xyz.h"
#include <chrono>
#include <cmath>
#include <iostream>
#include <thermostat.h>

const double DT{0.001};
const double SIGMA{1.44};
const double CUTOFF{5. * SIGMA};
// temperature in kelvins
const double TEMPERATURE{100.};
// number of times each run is repeated to get statistics about the runtime
// (must be >1 for Bessel's correction)
const size_t NUMBER_OF_RUNS{10};
// number of timesteps in each simulation
const size_t NUMBER_OF_TIMESTEPS{50};
// maximum number of atoms in the test series
const size_t NB_ATOMS_MAX{1000};
// step size of number of atoms in the test series
const size_t NB_ATOMS_STEP{50};

/// @brief Measure the execution time for the equilibration of a Lennard Jones
/// lattice with direct summation for a given number of atoms and simulation
/// steps
/// @param nb_atoms number of atoms in the lattice
/// @param direct whether to use direct summation or the truncated and shifted
/// potential
/// @return number of microseconds of execution time
double run_timed(size_t nb_atoms, bool direct) {
    Atoms atoms{Atoms(nb_atoms, SIGMA * pow(2, 1. / 6.))};
    // time the execution of a simulation as seen in milestone 5
    auto start = std::chrono::high_resolution_clock::now();
    (void)(direct ? lj_direct_summation(atoms, 1.0, SIGMA)
                  : ljts(atoms, 1.0, SIGMA, CUTOFF));
    for (size_t i{0}; i < NUMBER_OF_TIMESTEPS; i++) {
        verlet_step1(atoms.positions, atoms.velocities, atoms.forces, DT,
                     atoms.masses);
        (void)(direct ? lj_direct_summation(atoms, 1.0, SIGMA)
                      : ljts(atoms, 1.0, SIGMA, CUTOFF));
        verlet_step2(atoms.velocities, atoms.forces, DT, atoms.masses);
        berendsen_thermostat(atoms, TEMPERATURE, DT, 100 * DT);
    }
    auto stop = std::chrono::high_resolution_clock::now();
    // compute the duration of the loop
    return (double)std::chrono::duration_cast<std::chrono::microseconds>(stop -
                                                                         start)
        .count();
}

int main(int argc, char *argv[]) {
    // open a csv file and write the header describing the stored data
    std::ofstream file("runtimes.csv");
    // output all relevant information to the csv file in the format:
    // direct?,nb_atoms,average,min,max,stddev,runtime1, runtime2,...
    file << "direct summation or ljts,number of atoms,average runtime,minimum "
            "runtime,maximum runtime,corrected sample standard deviation";
    for (size_t i{1}; i <= NUMBER_OF_RUNS; i++) {
        file << ",runtime for run nr. " << i;
    }
    file << std::endl;

    // for direct summation and ljts and each desired number of atoms, run
    // NUMBER_OF_RUNS many simulations and collect statistical data on the
    // runtimes
    for (bool direct : {false, true}) {
        for (size_t nb_atoms{std::max((size_t)2, NB_ATOMS_STEP)};
             nb_atoms <= NB_ATOMS_MAX; nb_atoms += NB_ATOMS_STEP) {
            // output to stddout to keep track of the process
            std::cout << "Running with " << nb_atoms << " atoms" << std::endl;
            // the core loop where simulation runs are timed:
            double runs[NUMBER_OF_RUNS];
            for (size_t run = 0; run < NUMBER_OF_RUNS; run++) {
                runs[run] = run_timed(nb_atoms, direct);
            }
            // compute minimum, maximum and average runtime
            double average{0.};
            double min{__DBL_MAX__};
            double max{0.};
            for (double val : runs) {
                average += val;
                min = val < min ? val : min;
                max = val > max ? val : max;
            }
            average /= NUMBER_OF_RUNS;
            // compute the standard deviation with Bessel's correction
            double stddev{0.};
            for (double val : runs) {
                stddev += (val - average) * (val - average);
            }
            stddev = sqrt(stddev / (NUMBER_OF_RUNS - 1.));
            // write to csv in the same format as specified in the header:
            // direct,nb_atoms,average,min,max,stddev,runtime1, runtime2,...
            file << (direct ? "direct," : "ljts,") << nb_atoms << "," << average
                 << "," << min << "," << max << "," << stddev;
            for (size_t i{0}; i < NUMBER_OF_RUNS; i++) {
                file << "," << runs[i];
            }
            file << std::endl;
        }
    }
    return 0;
}