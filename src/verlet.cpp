#include "verlet.h"

/// The predictor step of a Velocity-Verlet time integration scheme
void verlet_step1(double &x, double &y, double &z, double &vx, double &vy,
                  double &vz, double fx, double fy, double fz, double dt,
                  double m) {
    // update each velocity component by a half step to obtain v(t + 0.5dt)
    vx += fx / m * dt * .5;
    vy += fy / m * dt * .5;
    vz += fz / m * dt * .5;
    // use v(t + 0.5dt) to perform a full explicit Euler step, obtaining a new
    // position
    x += vx * dt;
    y += vy * dt;
    z += vz * dt;
}

/// The corrector step of a velocity-Verlet time integration scheme.
/// Use this after verlet_step1 AND a subsequent update to fx,fy,fz using
/// the positions obtained through verlet_step1
void verlet_step2(double &vx, double &vy, double &vz, double fx, double fy,
                  double fz, double dt, double m) {
    // obtain the velocities v(t+dt) using the newly calculated forces at
    // r(t+dt) this is like another half-step of explicit euler but at the new
    // locations
    vx += fx / m * dt * .5;
    vy += fy / m * dt * .5;
    vz += fz / m * dt * .5;
}
