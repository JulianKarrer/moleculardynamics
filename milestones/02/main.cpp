#include "verlet.h"
#include <cmath>
#include <iostream>

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

/// A test scenario for velocity verlet integration where a point-mass orbits a
/// massive object at the origin. The initial velocity and position are chosen
/// such that a perfectly circular orbit in the xy-plane should be realised.
/// Deviations from this orbit are printed out for each timestep.
int main(int argc, char *argv[]) {
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
    double dt = 1.;
    // iteratively check assertions
    for (uint i = 0; i < 500; i++) {
        // calculate forces and integrate time
        gravity_forces_orbit_xy(x, y, z, fx, fy, fz, g, m);
        verlet_step1(x, y, z, vx, vy, vz, fx, fy, fz, dt, m);
        gravity_forces_orbit_xy(x, y, z, fx, fy, fz, g, m);
        verlet_step2(vx, vy, vz, fx, fy, fz, dt, m);
        std::cout << "Timestep " << i << "\t Radius deviation:"
                  << abs(r_0 - sqrt(x * x + y * y + z * z))
                  << "\t Position: x=" << x << " y=" << y << "\n";
    }
    std::cout << "(Deviation should be ~zero since energy is conserved)\n";
    return 0;
}
