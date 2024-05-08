#ifndef __ATOMS_H
#define __ATOMS_H

#include "types.h"

class Atoms {
  public:
    Positions_t positions;
    Velocities_t velocities;
    Forces_t forces;

    Atoms(const Positions_t &p);

    Atoms(const Positions_t &p, const Velocities_t &v);

    size_t nb_atoms() const;
};

#endif // __ATOMS_H