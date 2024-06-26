#ifndef __ATOMS_H
#define __ATOMS_H

#include "types.h"

/// @brief  A container class for atoms, holding relevant information such as
/// positions, velocities, currently acting forces, masses etc. Provides
/// convenience functions for querying properties of the system such as kinetic
/// energy or atom count.
class Atoms {
  public:
    Positions_t positions;
    Velocities_t velocities;
    Forces_t forces;
    Masses_t masses;
    Names_t names;

    Atoms(const size_t n);

    Atoms(const size_t n, double radius);

    Atoms(const Names_t names, Positions_t &p);

    Atoms(const Positions_t &p);

    Atoms(const Positions_t &p, const Velocities_t &v);

    Atoms(const Positions_t &p, const Velocities_t &v, const Masses_t &m);

    Atoms(const Positions_t &p, const Velocities_t &v, const Names_t names);

    size_t nb_atoms() const;

    double kinetic_energy();
};

/// @brief As described in
/// https://pastewka.github.io/MolecularDynamics/_notes/eV_A_units, the time
/// units are fixed at fs, energies at eV and lengths at Angstrom, such that
/// masses are given in units of 9.64853321...g/mol to make all units
/// consistent.
const double U_PER_MASS_UNIT{1.0 / (9.64853321e-6)};

/// Map names of elements to their repective atomic masses in units of u.
/// Data from https://pubchem.ncbi.nlm.nih.gov/periodic-table/
const std::unordered_map<std::string, double> ELEM_NAME_TO_MASS{
    {"H", U_PER_MASS_UNIT * 1.0080},       {"He", U_PER_MASS_UNIT * 4.00260},
    {"Li", U_PER_MASS_UNIT * 7.0},         {"Be", U_PER_MASS_UNIT * 9.012183},
    {"B", U_PER_MASS_UNIT * 10.81},        {"C", U_PER_MASS_UNIT * 12.011},
    {"N", U_PER_MASS_UNIT * 14.007},       {"O", U_PER_MASS_UNIT * 15.999},
    {"F", U_PER_MASS_UNIT * 18.99840316},  {"Ne", U_PER_MASS_UNIT * 20.180},
    {"Na", U_PER_MASS_UNIT * 22.9897693},  {"Mg", U_PER_MASS_UNIT * 24.305},
    {"Al", U_PER_MASS_UNIT * 26.981538},   {"Si", U_PER_MASS_UNIT * 28.085},
    {"P", U_PER_MASS_UNIT * 30.97376200},  {"S", U_PER_MASS_UNIT * 32.07},
    {"Cl", U_PER_MASS_UNIT * 35.45},       {"Ar", U_PER_MASS_UNIT * 39.9},
    {"K", U_PER_MASS_UNIT * 39.0983},      {"Ca", U_PER_MASS_UNIT * 40.08},
    {"Sc", U_PER_MASS_UNIT * 44.95591},    {"Ti", U_PER_MASS_UNIT * 47.867},
    {"V", U_PER_MASS_UNIT * 50.9415},      {"Cr", U_PER_MASS_UNIT * 51.996},
    {"Mn", U_PER_MASS_UNIT * 54.93804},    {"Fe", U_PER_MASS_UNIT * 55.84},
    {"Co", U_PER_MASS_UNIT * 58.93319},    {"Ni", U_PER_MASS_UNIT * 58.693},
    {"Cu", U_PER_MASS_UNIT * 63.55},       {"Zn", U_PER_MASS_UNIT * 65.4},
    {"Ga", U_PER_MASS_UNIT * 69.723},      {"Ge", U_PER_MASS_UNIT * 72.63},
    {"As", U_PER_MASS_UNIT * 74.92159},    {"Se", U_PER_MASS_UNIT * 78.97},
    {"Br", U_PER_MASS_UNIT * 79.90},       {"Kr", U_PER_MASS_UNIT * 83.80},
    {"Rb", U_PER_MASS_UNIT * 85.468},      {"Sr", U_PER_MASS_UNIT * 87.62},
    {"Y", U_PER_MASS_UNIT * 88.90584},     {"Zr", U_PER_MASS_UNIT * 91.22},
    {"Nb", U_PER_MASS_UNIT * 92.90637},    {"Mo", U_PER_MASS_UNIT * 95.95},
    {"Tc", U_PER_MASS_UNIT * 96.90636},    {"Ru", U_PER_MASS_UNIT * 101.1},
    {"Rh", U_PER_MASS_UNIT * 102.9055},    {"Pd", U_PER_MASS_UNIT * 106.42},
    {"Ag", U_PER_MASS_UNIT * 107.868},     {"Cd", U_PER_MASS_UNIT * 112.41},
    {"In", U_PER_MASS_UNIT * 114.818},     {"Sn", U_PER_MASS_UNIT * 118.71},
    {"Sb", U_PER_MASS_UNIT * 121.760},     {"Te", U_PER_MASS_UNIT * 127.6},
    {"I", U_PER_MASS_UNIT * 126.9045},     {"Xe", U_PER_MASS_UNIT * 131.29},
    {"Cs", U_PER_MASS_UNIT * 132.9054520}, {"Ba", U_PER_MASS_UNIT * 137.33},
    {"La", U_PER_MASS_UNIT * 138.9055},    {"Ce", U_PER_MASS_UNIT * 140.116},
    {"Pr", U_PER_MASS_UNIT * 140.90766},   {"Nd", U_PER_MASS_UNIT * 144.24},
    {"Pm", U_PER_MASS_UNIT * 144.91276},   {"Sm", U_PER_MASS_UNIT * 150.4},
    {"Eu", U_PER_MASS_UNIT * 151.964},     {"Gd", U_PER_MASS_UNIT * 157.2},
    {"Tb", U_PER_MASS_UNIT * 158.92535},   {"Dy", U_PER_MASS_UNIT * 162.500},
    {"Ho", U_PER_MASS_UNIT * 164.93033},   {"Er", U_PER_MASS_UNIT * 167.26},
    {"Tm", U_PER_MASS_UNIT * 168.93422},   {"Yb", U_PER_MASS_UNIT * 173.05},
    {"Lu", U_PER_MASS_UNIT * 174.9668},    {"Hf", U_PER_MASS_UNIT * 178.49},
    {"Ta", U_PER_MASS_UNIT * 180.9479},    {"W", U_PER_MASS_UNIT * 183.84},
    {"Re", U_PER_MASS_UNIT * 186.207},     {"Os", U_PER_MASS_UNIT * 190.2},
    {"Ir", U_PER_MASS_UNIT * 192.22},      {"Pt", U_PER_MASS_UNIT * 195.08},
    {"Au", U_PER_MASS_UNIT * 196.96657},   {"Hg", U_PER_MASS_UNIT * 200.59},
    {"Tl", U_PER_MASS_UNIT * 204.383},     {"Pb", U_PER_MASS_UNIT * 207},
    {"Bi", U_PER_MASS_UNIT * 208.98040},   {"Po", U_PER_MASS_UNIT * 208.98243},
    {"At", U_PER_MASS_UNIT * 209.98715},   {"Rn", U_PER_MASS_UNIT * 222.01758},
    {"Fr", U_PER_MASS_UNIT * 223.01973},   {"Ra", U_PER_MASS_UNIT * 226.02541},
    {"Ac", U_PER_MASS_UNIT * 227.02775},   {"Th", U_PER_MASS_UNIT * 232.038},
    {"Pa", U_PER_MASS_UNIT * 231.03588},   {"U", U_PER_MASS_UNIT * 238.0289},
    {"Np", U_PER_MASS_UNIT * 237.048172},  {"Pu", U_PER_MASS_UNIT * 244.06420},
    {"Am", U_PER_MASS_UNIT * 243.061380},  {"Cm", U_PER_MASS_UNIT * 247.07035},
    {"Bk", U_PER_MASS_UNIT * 247.07031},   {"Cf", U_PER_MASS_UNIT * 251.07959},
    {"Es", U_PER_MASS_UNIT * 252.0830},    {"Fm", U_PER_MASS_UNIT * 257.09511},
    {"Md", U_PER_MASS_UNIT * 258.09843},   {"No", U_PER_MASS_UNIT * 259.10100},
    {"Lr", U_PER_MASS_UNIT * 266.120},     {"Rf", U_PER_MASS_UNIT * 267.122},
    {"Db", U_PER_MASS_UNIT * 268.126},     {"Sg", U_PER_MASS_UNIT * 269.128},
    {"Bh", U_PER_MASS_UNIT * 270.133},     {"Hs", U_PER_MASS_UNIT * 269.1336},
    {"Mt", U_PER_MASS_UNIT * 277.154},     {"Ds", U_PER_MASS_UNIT * 282.166},
    {"Rg", U_PER_MASS_UNIT * 282.169},     {"Cn", U_PER_MASS_UNIT * 286.179},
    {"Nh", U_PER_MASS_UNIT * 286.182},     {"Fl", U_PER_MASS_UNIT * 290.192},
    {"Mc", U_PER_MASS_UNIT * 290.196},     {"Lv", U_PER_MASS_UNIT * 293.205},
    {"Ts", U_PER_MASS_UNIT * 294.211},     {"Og", U_PER_MASS_UNIT * 295.216},
};

#endif // __ATOMS_H