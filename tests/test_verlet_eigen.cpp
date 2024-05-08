#include "verlet.h"
#include <Eigen/Dense>
#include <cmath>
#include <gtest/gtest.h>
#include <iostream>

// Test if all point masses under constant force of gravity
// 1) accelerate with 9.81m/s^2, ie. |vy| increases by dt*g per timestep
// 2) follow the analytic solution s = 1/2 * g * t^2
TEST(VerletEigenTest, GravityParticleAcceleratesCorrectly) {
    // intialize state of the system
    int nb_atoms = 10;
    Positions_t positions(3, nb_atoms);
    Velocities_t velocities(3, nb_atoms);
    Velocities_t prev_velocities(3, nb_atoms);
    Forces_t forces(3, nb_atoms);
    double m = 1.;
    double dt = .1;
    double g = -9.81;
    double t = 0.;
    positions.setZero();
    velocities.setZero();
    prev_velocities.setZero();
    for (int i = 0; i < nb_atoms; i++) {
        forces.col(i) = Vec3_t{0., g * m, 0.};
    }
    for (uint step = 0; step < 100; step++) {
        verlet_step1(positions, velocities, forces, dt, m);
        verlet_step2(velocities, forces, dt, m);
        t += dt;
        for (int i = 0; i < nb_atoms; i++) {
            // assertion 1) vy decreases by 0.981/timestep
            ASSERT_NEAR(abs(velocities(1, i)),
                        abs(prev_velocities(1, i) + dt * g), 1e-8);
            // 2) assert equivalence to analytic solution s = 1/2 * g * t^2 for
            // each particle
            ASSERT_NEAR(positions(1, i), 0.5 * g * t * t, 1e-8);
        }

        prev_velocities = velocities;
    }
}

// Test if all point masses under gravity at the right initial speed
// (|v_init| = sqrt(g/r_0)) orbit a mass at the origin
// 1) in a perfect circle
// 2) at constant velocity
// 3) only in the initial plane of the orbit

inline void gravity_forces_orbit_xy(Positions_t positions, Forces_t &forces,
                                    double g, double m, int nb_atoms) {
    for (int i = 0; i < nb_atoms; i++) {
        double x = positions(0, i);
        double y = positions(1, i);
        double z = positions(2, i);
        double r = sqrt(x * x + y * y + z * z);
        forces.col(i) = Vec3_t{(x / r) * (-g * m / (r * r)) * m,
                               (y / r) * (-g * m / (r * r)) * m,
                               (z / r) * (-g * m / (r * r)) * m};
    }
}

TEST(VerletEigenTest, GravityParticleOrbits) {
    // intialize state of the system
    int nb_atoms = 100;
    Positions_t positions(3, nb_atoms);
    Velocities_t velocities(3, nb_atoms);
    Forces_t forces(3, nb_atoms);
    double m = 1.;
    double dt = 0.1;
    double g = 100;
    // spawn atoms at radii 150.+2.*i
    forces.setZero();
    for (int i = 0; i < nb_atoms; i++) {
        double r = 150. + 2. * ((double)i);
        positions.col(i) = Vec3_t{r, 0., 0.};
        velocities.col(i) = Vec3_t{0., -sqrt(g / r), 0.};
    }
    for (uint step = 0; step < 500; step++) {
        gravity_forces_orbit_xy(positions, forces, g, m, nb_atoms);
        verlet_step1(positions, velocities, forces, dt, m);
        gravity_forces_orbit_xy(positions, forces, g, m, nb_atoms);
        verlet_step2(velocities, forces, dt, m);
        for (int i = 0; i < nb_atoms; i++) {
            double r = 150. + 2. * ((double)i);
            double x = positions(0, i);
            double y = positions(1, i);
            double z = positions(2, i);
            ASSERT_NEAR(sqrt(x * x + y * y + z * z), r, 1e-3);
        }
    }
}