#ifndef __VERLET_H
#define __VERLET_H

#include <Eigen/Dense>
#include <types.h>

void verlet_step1(double &x, double &y, double &z, double &vx, double &vy,
                  double &vz, double fx, double fy, double fz, double dt,
                  double m);
void verlet_step2(double &vx, double &vy, double &vz, double fx, double fy,
                  double fz, double dt, double m);
void verlet_step1(Positions_t &positions, Velocities_t &velocities,
                  const Forces_t &forces, double dt, Masses_t m);
void verlet_step2(Velocities_t &velocities, const Forces_t &forces, double dt,
                  Masses_t m);
void verlet_step1(Positions_t &positions, Velocities_t &velocities,
                  const Forces_t &forces, double dt, double m);
void verlet_step2(Velocities_t &velocities, const Forces_t &forces, double dt,
                  double m);
#endif // __VERLET_H
