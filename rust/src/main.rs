use std::{
    f64::consts::SQRT_2,
    fs::{File, OpenOptions},
    io::{BufWriter, Write},
    time::Instant,
};

use kiddo::{
    immutable::float::kdtree::ImmutableKdTree, KdTree, NearestNeighbour, SquaredEuclidean,
};
use nalgebra::Vector3;
use rayon::iter::{
    IndexedParallelIterator, IntoParallelRefIterator, IntoParallelRefMutIterator, ParallelIterator,
};

// switch default allocator
use mimalloc::MiMalloc;
#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;

// Define types
/// Type alias for typechecking positions
type PositionT = Vec<Vector3<f64>>;
/// Type alias for typechecking velocities
type VelocityT = Vec<Vector3<f64>>;
/// Type alias for typechecking forces
type ForceT = Vec<Vector3<f64>>;
/// Type alias for typechecking masses
type MassT = Vec<f64>;

// Define constants
/// The Boltzmann constant defined in terms of Kelvin and electron Volts
const KB_EV: f64 = 8.617333262e-5;

/// Structure that holds buffers for the positions, velocities, forces, masses etc. of all atoms.
/// We use a structure of arrays instead of an array of structs for performance reasons.
struct Atoms {
    pos: PositionT,
    vel: VelocityT,
    frc: ForceT,
    mas: MassT,
}

impl Atoms {
    /// Return the number of atoms as an f64 for convenience
    pub fn nb_atoms(&self) -> f64 {
        self.pos.len() as f64
    }

    /// Initialize a set of n grid points with given spacing on an oblique lattice
    fn initialize_lattice(nb_atoms: usize, spacing: f64) -> PositionT {
        let cube_length: usize = (nb_atoms as f64).cbrt().ceil() as usize;
        // offset the whole structure by half the estimated size to centre around the origin
        let centre = cube_length as f64 * spacing * 0.5;
        let mut res: PositionT = vec![Vector3::new(0.0, 0.0, 0.0); nb_atoms];
        let mut i = 0;
        for x in 0..cube_length + 1 {
            for y in 0..cube_length + 1 {
                for z in 0..cube_length + 1 {
                    if i >= nb_atoms {
                        break;
                    }
                    // offset every other xy plane by 0.5*spacing for an oblique lattice
                    let offset = if z % 2 == 0 {
                        0.25 * spacing
                    } else {
                        -0.25 * spacing
                    };
                    let d = i as f64 * 1e-15; // minimal offset so kdtree can fit planes
                    res[i] = Vector3::new(
                        x as f64 * spacing - centre + offset + d,
                        y as f64 * spacing - centre + offset + d,
                        z as f64 * spacing / SQRT_2 - centre / SQRT_2 + d,
                    );
                    i += 1;
                }
            }
        }
        res
    }

    /// Initialize a set of `nb_atoms` on an oblique lattice with unit mass and no velocities or forces
    pub fn new(nb_atoms: usize, sigma: f64) -> Self {
        let pos = Atoms::initialize_lattice(nb_atoms, sigma * 2.0f64.powf(1. / 6.));
        return Self {
            pos,
            vel: vec![Vector3::zeros(); nb_atoms],
            frc: vec![Vector3::zeros(); nb_atoms],
            mas: vec![1.0; nb_atoms],
        };
    }

    pub fn lj_direct_sum(&mut self, epsilon: f64, sigma: f64) -> f64 {
        self.frc.fill(Vector3::zeros());
        self.pos
            .iter()
            .enumerate()
            .map(|(i, x_i)| {
                self.pos
                    .iter()
                    .enumerate()
                    .map(|(j, x_j)| {
                        if i < j {
                            let r = (x_i - x_j).norm(); // sqrt of the squared distance
                            let r_ij_hat: Vector3<f64> =
                                (x_i - x_j).try_normalize(1e-15).unwrap_or(Vector3::zeros());
                            let f_ij: Vector3<f64> = r_ij_hat * lj_pot_deriv(r, epsilon, sigma);
                            self.frc[i] -= f_ij;
                            self.frc[j] += f_ij;
                            // return the shifted and truncated potential
                            lj_pot(r, epsilon, sigma)
                        } else {
                            0.0
                        }
                    })
                    .sum::<f64>()
            })
            .sum::<f64>()
    }

    pub fn ljts(&mut self, epsilon: f64, sigma: f64, cutoff: f64) -> f64 {
        // reset forces, since we add later!
        self.frc.fill(Vector3::zeros());
        // create immutable kd tree, unfortunately requires taking ownership of the buffer/copying
        // TODO: check if it is possible to use slices/views instead
        let kd = kiddo::ImmutableKdTree::new_from_slice(
            &self
                .pos
                .iter()
                .map(|p| p.as_ref().to_owned())
                .collect::<Vec<[f64; 3]>>(),
        );
        // update forces and return the sum of all pair-potentials
        self.pos
            .iter()
            .enumerate()
            .map(|(i, x_i)| {
                kd.within::<SquaredEuclidean>(&[x_i.x, x_i.y, x_i.z], cutoff)
                    .iter()
                    .map(|NearestNeighbour { distance, item }| {
                        // compute force
                        let r = distance.sqrt(); // sqrt of the squared distance
                        let j = *item as usize;
                        if i < j && r > 0.0 {
                            let x_j = self.pos[j];
                            let r_ij_hat: Vector3<f64> =
                                (x_i - x_j).try_normalize(1e-15).unwrap_or(Vector3::zeros());
                            let f_ij: Vector3<f64> = r_ij_hat * lj_pot_deriv(r, epsilon, sigma);
                            self.frc[i] -= f_ij;
                            self.frc[j] += f_ij;
                            // return the shifted and truncated potential
                            ljts_pot(r, epsilon, sigma, cutoff)
                        } else {
                            0.0
                        }
                    })
                    .sum::<f64>()
            })
            .sum()
    }

    pub fn kinetic_energy(&self) -> f64 {
        self.vel
            .iter()
            .zip(&self.mas)
            .map(|(v, m)| v.norm_squared() * m * 0.5)
            .sum()
    }

    pub fn temperature(&self) -> f64 {
        self.kinetic_energy() / (KB_EV * 1.5 * self.nb_atoms())
    }

    pub fn berendsen_thermostat(&mut self, target_temp: f64, dt: f64, relaxation: f64) {
        let current_temp = self.temperature();
        // avoid division by zero
        if current_temp > 0.0 {
            self.vel.iter_mut().for_each(|v| {
                let lambda = (1. + (target_temp / current_temp - 1.) * dt / relaxation).sqrt();
                *v *= lambda;
            })
        }
    }

    /// Write current positions and velocities to .xyz file
    pub fn write_to_xyz(&self, file: &mut BufWriter<File>) {
        write!(file, "{}\n\n", self.nb_atoms()).unwrap();
        self.pos.iter().zip(&self.vel).for_each(|(x, v)| {
            write!(
                file,
                "H\t{}\t{}\t{}\t{}\t{}\t{}\t\n",
                x.x, x.y, x.z, v.x, v.y, v.z
            )
            .unwrap();
        });
    }
}

/// The Lennard-Jones potential
fn lj_pot(r: f64, epsilon: f64, sigma: f64) -> f64 {
    4. * epsilon * ((sigma / r).powi(12) - (sigma / r).powi(6))
}

/// The truncated and shifted Lennard-Jones potential with given cutoff distance
fn ljts_pot(r: f64, epsilon: f64, sigma: f64, cutoff: f64) -> f64 {
    if r > cutoff {
        0.0
    } else {
        lj_pot(r, epsilon, sigma) - lj_pot(cutoff, epsilon, sigma)
    }
}

/// Magnitude of the derivative of the Lennard-Jones potential
fn lj_pot_deriv(r: f64, epsilon: f64, sigma: f64) -> f64 {
    let sig_over_r = sigma / r;
    epsilon * (-48. * sig_over_r.powi(12) / r + 24. * sig_over_r.powi(6) / r)
}

/// Predictor step of the Velocity Verlet time integration scheme
fn velocity_verlet_step1(
    pos: &mut PositionT,
    vel: &mut VelocityT,
    frc: &ForceT,
    mas: &MassT,
    dt: f64,
) {
    vel.iter_mut()
        .zip(pos)
        .zip(frc)
        .zip(mas)
        .for_each(|(((v, x), f), m)| {
            *v += (f * (1. / *m)) * dt * 0.5;
            *x += *v * dt;
        });
}

/// Corrector step of the Velocity Verlet time integration scheme
fn velocity_verlet_step2(vel: &mut VelocityT, frc: &ForceT, mas: &MassT, dt: f64) {
    vel.iter_mut().zip(frc).zip(mas).for_each(|((v, f), m)| {
        *v += (*f * (1. / m)) * dt * 0.5;
    });
}

// define parameters
const DT: f64 = 0.001;
const SIGMA: f64 = 1.44;
const EPSILON: f64 = 1.;
const CUTOFF: f64 = 5. * SIGMA;
const TEMPERATURE: f64 = 100.;

const NB_RUNS: usize = 10;
const NB_TIMESTEPS: usize = 50;
const NB_ATOMS_MAX: usize = 750;
const NB_ATOMS_STEP: usize = 50;

fn run_timed(nb_atoms: usize, direct: bool) -> f64 {
    let mut atoms = Atoms::new(nb_atoms, SIGMA);

    // time the execution of a simulation as seen in milestone 5
    let start = Instant::now();
    if direct {
        atoms.lj_direct_sum(EPSILON, SIGMA)
    } else {
        atoms.ljts(EPSILON, SIGMA, CUTOFF)
    };
    for _ in 0..NB_TIMESTEPS {
        // compute lj forces and integrate with velocity verlet
        velocity_verlet_step1(&mut atoms.pos, &mut atoms.vel, &atoms.frc, &atoms.mas, DT);
        if direct {
            atoms.lj_direct_sum(EPSILON, SIGMA)
        } else {
            atoms.ljts(EPSILON, SIGMA, CUTOFF)
        };
        velocity_verlet_step2(&mut atoms.vel, &atoms.frc, &atoms.mas, DT);

        // equilibrate the system with a berendesen thermostat
        atoms.berendsen_thermostat(TEMPERATURE, DT, 100. * DT);
    }
    let time_taken = start.elapsed().as_micros();
    time_taken as f64
}

fn main() {
    // open trajectory file
    let mut file: BufWriter<File> = BufWriter::new(
        OpenOptions::new()
            .write(true)
            .create(true)
            .truncate(true)
            .open("runtimes.csv")
            .unwrap(),
    );
    // write csv header
    write!(file, "direct summation or ljts,number of atoms,average runtime,minimum runtime,maximum runtime,corrected sample standard deviation\n").unwrap();

    for direct in [false, true] {
        for n in (2.max(NB_ATOMS_STEP)..=NB_ATOMS_MAX).step_by(NB_ATOMS_STEP) {
            println!("Running with {} atoms", n);
            let runs: Vec<f64> = (0..NB_RUNS).map(|_| run_timed(n, direct)).collect();
            let average = runs.iter().sum::<f64>() / (NB_RUNS as f64);
            let min = runs.iter().min_by(|x, y| x.total_cmp(y)).unwrap();
            let max = runs.iter().max_by(|x, y| x.total_cmp(y)).unwrap();
            let stddev = (runs.iter().map(|val| (val - average).powi(2)).sum::<f64>()
                / (NB_RUNS as f64 - 1.))
                .sqrt();
            // write data to csv
            writeln!(
                file,
                "{},{},{},{},{},{}",
                if direct { "direct" } else { "ljts" },
                n,
                average,
                min,
                max,
                stddev
            )
            .unwrap();
        }
    }

    // main simulation loop
}
