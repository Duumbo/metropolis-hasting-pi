/// Compute pi with Metropolis-Hasting
use rand::prelude::*;
use std::fs::File;
use std::io::prelude::*;

// Grid size
const SIZE: usize = 10;
// Radius of the inscribed circle
const RADIUS: usize = SIZE / 2;
const NSAMP: usize = 10_000;
const NWARM: usize = 0;
// Sample interval
const NSKIP: usize = 10;
const N_RUNS: usize = 1;
// If true, probabilité d'acceptation = 0.5 * (n_voisins new) / (n_voisins old)
// If false, probabilité d'acceptation = 0.5
const NEIGHBORS_CHECK: bool = true;
const SEED: u64 = 101203;
const INITIAL_POINT: (usize, usize) = (0, 0);

fn compute_n_neighbors(coords: (usize, usize)) -> usize {
    let x1 = coords.0;
    let x2 = coords.1;
    let mut neighbors = 8;
    if x1 == 0 {
        neighbors -= 3;
    }
    if x2 == 0 {
        neighbors -= 3;
    }
    if x1 == SIZE - 1 {
        neighbors -= 3;
    }
    if x2 == SIZE - 1 {
        neighbors -= 3;
    }
    let boundary_cond: bool = (x1 == 0 && x2 == 0) ||
        (x1 == 0 && x2 == SIZE - 1) ||
        (x1 == SIZE - 1 && x2 == 0) ||
        (x1 == SIZE - 1 && x2 == SIZE - 1);
    if boundary_cond {
        neighbors = 3;
    }
    neighbors
}

fn is_legal(x1: usize, x2: usize, rngx: usize, rngy: usize) -> bool {
    if x1 == 0 {
        if rngx == 0 {
            return false;
        }
    }
    if x2 == 0 {
        if rngy == 0 {
            return false;
        }
    }
    if x1 == SIZE - 1 {
        if rngx == 2 {
            return false;
        }
    }
    if x2 == SIZE - 1 {
        if rngy == 2 {
            return false;
        }
    }
    if rngx == 1 && rngy == 1 {
        return false;
    }
    return true;
}

fn propose<R>(x: (usize, usize), rng: &mut R) -> (usize, usize)
    where R: Rng
{
    let mut x1 = x.0;
    let mut x2 = x.1;
    let mut coords_out = x;
    let mut legal_guess = false;
    let previous_neighbors = compute_n_neighbors(x);
    while ! legal_guess {
        let rngx = (rng.random::<u32>() % 3) as usize;
        let rngy = (rng.random::<u32>() % 3) as usize;
        if ! is_legal(x1, x2, rngx, rngy) {
            continue;
        }

        // Update coords
        x1 += rngx;
        x1 -= 1;
        x2 += rngy;
        x2 -= 1;

        let new_neighbors = compute_n_neighbors((x1, x2));
        let neighbors_ratio = previous_neighbors as f64 / new_neighbors as f64;
        if (rng.random::<f64>() < neighbors_ratio) || !NEIGHBORS_CHECK {
            // Propose
            coords_out = (x1, x2);
        }
        legal_guess = true;
    }
    coords_out
}

fn compute_pi() -> f64 {
    // Our sampling point
    let mut x = INITIAL_POINT;
    let mut count = 0;
    let mut fp = File::create("distribution").unwrap();
    fp.write(format!("# {} {} {} ", SIZE, NSAMP, NWARM).as_bytes()).unwrap();

    let mut rng = SmallRng::seed_from_u64(SEED);
    for _ in 0..NWARM {
        let xprop = propose(x, &mut rng);
        if rng.random_bool(0.5) {
            // Accept
            x = xprop;
        }
    }

    println!("Initial point = ({}, {})", x.0, x.1);
    for i in 0..NSAMP * NSKIP {
        let xprop = propose(x, &mut rng);
        if rng.random_bool(0.5) {
            // Accept
            x = xprop;
        }

        // Distance from origin
        let (xrel, yrel) = (x.0 as isize - RADIUS as isize, x.1 as isize - RADIUS as isize);
        let distsq = (xrel * xrel + yrel * yrel) as usize;

        if i % NSKIP == 0 {
            // Sample
            fp.write(format!("{}\n", x.0 + SIZE * x.1).as_bytes()).unwrap();
            if distsq <= RADIUS * RADIUS {
                count += 1;
            }
        }
    }
    let _area_sq = NSAMP;
    let _area_circ = count;
    // A_s / A_c = pi r^2 / 4 r^2
    println!("count = {}", count);
    let pi = 4.0 * count as f64 / (NSAMP) as f64;
    println!("pi = {}", pi);
    pi
}

fn main() {
    let mut pi_avg = 0.0;
    let mut pi_dev = 0.0;
    for _ in 0..N_RUNS {
        let pi = compute_pi();
        pi_avg += pi;
        pi_dev += pi*pi;
    }
    pi_avg = pi_avg / N_RUNS as f64;
    pi_dev = pi_dev / N_RUNS as f64 - pi_avg*pi_avg;
    println!("Avg = {}", pi_avg);
    println!("Dev = {}", pi_dev);
}
