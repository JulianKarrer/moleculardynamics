#ifndef __VERLET_H
#define __VERLET_H

void verlet_step1(double &x, double &y, double &z, double &vx, double &vy, double &vz,
                  double fx, double fy, double fz, double dt, double m);
void verlet_step2(double &vx, double &vy, double &vz, double fx, double fy, double fz,
                  double dt, double m);

#endif  // __VERLET_H
