// Include a library file to make sure proper includes are set
#include "verlet.h"
#include <cmath>
#include <gtest/gtest.h>
#include <iostream>

// Test if a point mass under constant force of gravity
// 1) does in fact fall downwards
// 2) with no sidewards motion
// 3) and increasing speed in y direction
// 4) but no other direction
TEST(VerletTest, GravityParticleFalls) {
    // intialize state of the system
    double px, py, pz, pvx, pvy, pvz, x, y, z, vx, vy, vz, fx, fz;
    px = py = pz = pvx = pvy = pvz = x = y = z = vx = vy = vz = fx = fz = 0.;
    double fy = -9.81;
    double m = 1.;
    double dt = 0.01;
    // iteratively check assertions
    for (uint i = 0; i < 50; i++) {
        verlet_step1(x, y, z, vx, vy, vz, fx, fy, fz, dt, m);
        verlet_step2(vx, vy, vz, fx, fy, fz, dt, m);
        // assertion 1) - falling down
        ASSERT_LT(y, py);
        // assertion 2) - not sidewards
        ASSERT_NEAR(x, px, 1e-8);
        ASSERT_NEAR(z, pz, 1e-8);
        // assertion 3) - increasing speed downwards
        ASSERT_LE(abs(pvy), abs(vy));
        ASSERT_LE(py, 0.);
        // assertion 4) - not sidewards
        ASSERT_NEAR(vx, pvx, 1e-8);
        ASSERT_NEAR(vz, pvz, 1e-8);
        // px,py,pz,pvx,pvy,pvz hold previous positions and velocities
        // from the past loop iteration
        px = x;
        py = y;
        pz = z;
        pvx = vx;
        pvy = vy;
        pvz = vz;
    }
}

// Test if a point mass under constant force of gravity
// 1) accelerates with 9.81m/s^2, ie. |vy| increases by dt*g per timestep
// 2) follows the analytic solution s = 1/2 * g * t^2
TEST(VerletTest, GravityParticleAcceleratesCorrectly) {
    // intialize state of the system
    double pvy, x, y, z, vx, vy, vz, fx, fz, t;
    pvy = x = y = z = vx = vy = vz = fx = fz = t = 0.;
    double fy, g;
    g = -9.81;
    double m = 1.;
    fy = g * m;
    double dt = 0.1;
    // iteratively check assertions
    for (uint i = 0; i < 50; i++) {
        verlet_step1(x, y, z, vx, vy, vz, fx, fy, fz, dt, m);
        verlet_step2(vx, vy, vz, fx, fy, fz, dt, m);
        // assertion 1) vy decreases by 0.981/timestep
        ASSERT_NEAR(abs(vy), abs(pvy + dt * g), 1e-8);
        // assertion 2) equivalence to analytic solution
        t += dt;
        ASSERT_NEAR(y, 0.5 * g * t * t, 1e-8);
        // store y-velocity from the previous loop iteration
        pvy = vy;
    }
}

/// @brief Define the force of gravity computed towards an unmoving mass at the
/// origin
/// @param x x component of the position
/// @param y y component of the position
/// @param z z component of the position
/// @param fx x component of the force computed
/// @param fy y component of the force computed
/// @param fz z component of the force computed
/// @param g the gravitational constant used
/// @param m the mass of the orbiting body
inline void gravity_forces_orbit_xy(double x, double y, double z, double &fx,
                                    double &fy, double &fz, double g,
                                    double m) {
    double r = sqrt(x * x + y * y + z * z);
    fx = (x / r) * (-g * m / (r * r)) * m;
    fy = (y / r) * (-g * m / (r * r)) * m;
    fz = (z / r) * (-g * m / (r * r)) * m;
}

// Test if a point mass under gravity at the right initial speed
// (|v_init| = sqrt(g/r_0)) orbits a mass at the origin
// 1) in a perfect circle
// 2) at constant velocity
// 3) only in the initial plane of the orbit
TEST(VerletTest, GravityParticleOrbits) {
    // intialize state of the system
    // gravitational constant
    double g = 100;
    // position
    double x, y, z;
    y = z = 0.;
    // radius
    double r_0 = 150.;
    x = r_0;
    // velocity components
    double vx, vy, vz;
    vx = vz = 0.;
    // initial velocity
    vy = -sqrt(g / r_0);
    // forces
    double fx, fy, fz;
    fx = fy = fz = 0.;
    // mass
    double m = 1.;
    // timestep
    double dt = 0.1;
    // iteratively check assertions
    gravity_forces_orbit_xy(x, y, z, fx, fy, fz, g, m);
    for (uint i = 0; i < 500; i++) {
        // calculate forces and integrate time
        verlet_step1(x, y, z, vx, vy, vz, fx, fy, fz, dt, m);
        gravity_forces_orbit_xy(x, y, z, fx, fy, fz, g, m);
        verlet_step2(vx, vy, vz, fx, fy, fz, dt, m);
        // 1)  in a perfect circle
        ASSERT_NEAR(sqrt(x * x + y * y + z * z), r_0, 1e-6);
        // 2) at constant velocity
        ASSERT_NEAR(sqrt(vx * vx + vy * vy + vz * vz), sqrt(g / r_0), 1e-8);
        // 3) only in the initial plane of the orbit
        ASSERT_NEAR(z, 0., 1e-8);
        ASSERT_NEAR(vz, 0., 1e-8);
    }
}