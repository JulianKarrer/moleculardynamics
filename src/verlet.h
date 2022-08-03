#ifndef VERLET_H
#define VERLET_H


void scalar_verlet_step1(double &x, double &y, double &z, double &vx,
                         double &vy, double &vz, double fx, double fy,
                         double fz, double timestep);
void scalar_verlet_step2(double &vx, double &vy, double &vz, double fx,
                         double fy, double fz, double timestep);

void nof();
#endif // __VERLET_H
