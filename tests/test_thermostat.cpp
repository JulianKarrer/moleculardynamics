#include <gtest/gtest.h>
#include <lj_direct_summation.h>
#include <thermostat.h>
#include <verlet.h>

/// @brief Implements the test of the thermostat for the different initial
/// conditions defined in the `TEST`s in the following
void berendsen_test(Atoms &atoms, double t_init, double t_wish, double dt,
                    double assert_tolerance) {
    // Perform 20k simulation steps to make sure the system should have
    // calibrated
    lj_direct_summation(atoms, 1., 1.);
    for (uint step = 0; step < 10000; step++) {
        verlet_step1(atoms.positions, atoms.velocities, atoms.forces, dt,
                     atoms.masses);
        lj_direct_summation(atoms, 1., 1.);
        verlet_step2(atoms.velocities, atoms.forces, dt, atoms.masses);
        berendsen_thermostat(atoms, t_wish, dt, 100 * dt);
    }
    double t_end{temperature_cur(atoms)};

    // Print results for debugging and to catch errors
    std::cout << "t_init: " << t_init << "t_wish: " << t_wish
              << "t_end: " << t_end << "\n";

    // Assert that the system has equilibrated.
    ASSERT_NEAR(t_end, t_wish, assert_tolerance);
}

/// @brief Tests if the Berendsen thermostat equlibrates a system of 10 atoms on
/// a grid
TEST(BerendsenThermosTest, Equilibrate) {
    // intialize state of the system
    Atoms atoms{Atoms(10, pow(2, 1. / 6.))};

    // Set target temperature and time step size
    double t_init{temperature_cur(atoms)};
    double t_wish{20.};
    double dt{0.001};

    berendsen_test(atoms, t_init, t_wish, dt, 0.1 * abs(t_wish - t_init));
}

/// @brief Check if the thermostat can cool the system to zero K. This should
/// result in standstill. The initially densly packed atoms should also expand,
/// meaning the bounding volume increases.
TEST(BerendsenThermosTest, ToZero) {
    Atoms atoms{Atoms(10, pow(2, 1. / 6.))};
    double initial_volume{(atoms.positions.rowwise().maxCoeff() -
                           atoms.positions.rowwise().minCoeff())
                              .prod()};
    double t_init{temperature_cur(atoms)};
    // set the target temperature to 0, otherwise perform the same test
    double t_wish{0};
    double dt{0.001};

    // very low tolerance, since this state should be easily reached
    berendsen_test(atoms, t_init, t_wish, dt, 1e-9);
    double final_volume{(atoms.positions.rowwise().maxCoeff() -
                         atoms.positions.rowwise().minCoeff())
                            .prod()};
    // assert that the bounding volume has expanded.
    std::cout << "initial volume: " << initial_volume
              << "final volume: " << final_volume << "\n";
    ASSERT_GT(final_volume, initial_volume);
}

/// @brief Check if a single atom at the origin with initial unit velocity
/// in x-direction speeds up to the velocity dictated by the target 
/// kinetic energy
TEST(BerendsenThermosTest, SingleAtom) {
    // set initial position at origin and initial velocity to unit vector in x-direction
    Atoms atoms{Atoms(1)};
    atoms.positions.setZero();
    atoms.velocities.col(0) = Vec3_t{1.,0.,0.};
    // set target 'temperature'/kinetic energy
    double t_init{temperature_cur(atoms)};
    double t_wish{t_init * 5.};
    double dt{0.001};

    // run the simulation
    berendsen_test(atoms, t_init, t_wish, dt, 1e-9);
    // assert the final kinetic energy is 5 times the initial kinetic energy:
    // should be initially 0.5*1*1*1=0.5, finally 5*0.5=2.5
    ASSERT_NEAR(atoms.kinetic_energy(), 2.5, 1e-9);
}