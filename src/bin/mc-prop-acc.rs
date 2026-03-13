/// Compute pi with Metropolis-Hasting
use rand::prelude::*;
use std::fs::File;
use std::io::prelude::*;
use std::sync::Mutex;
use std::thread;
use std::sync::Arc;

// Grid size
const SIZE: usize = 100;
// Radius of the inscribed circle
const RADIUS: usize = SIZE / 2;
const NSAMP: usize = 200;
const NWARM: usize = 1000;
// Sample interval
const NSKIP: usize = 1;
const N_RUNS: usize = 1000;
// If true, probabilité d'acceptation = 0.5 * (n_voisins new) / (n_voisins old)
// If false, probabilité d'acceptation = 0.5
const NEIGHBORS_CHECK: bool = true;
const SEED: u64 = 101203;
const INITIAL_POINT: (usize, usize) = (0, 0);
// Gaussian sigma and mean
const SIGMA: f64 = SIZE as f64 / 10.0;
const XMEAN: usize = SIZE/2;
const YMEAN: usize = SIZE/2;

fn rho(point: (usize, usize)) -> f64 {
    let shiftx = 1;
    let shifty = 1;
    let x = point.0 as f64 - (XMEAN -shiftx) as f64;
    let y = point.1 as f64 - (YMEAN-shifty) as f64;
    let c = 0.0;
    <f64>::exp( - (x*x + y*y) / (2.0 * (SIGMA+c)))
}

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

fn compute_pi(fp_mutex: Arc<Mutex<File>>, rng: u64) -> f64 {
    // Our sampling point
    let mut x = INITIAL_POINT;
    let mut count = 0;
    let mut sample = 0.0;

    let mut rng = SmallRng::seed_from_u64(SEED + rng);
    for _ in 0..NWARM {
        let xprop = propose(x, &mut rng);

        let rho1 = rho((x.0, x.1));
        let rho2 = rho((xprop.0, xprop.1));
        let ratio = rho2 / rho1;

        if rng.random_bool(ratio.min(1.0)) {
            // Accept
            x = xprop;
        }
    }

    println!("Initial point = ({}, {})", x.0, x.1);
    for i in 0..NSAMP * NSKIP {
        let xprop = propose(x, &mut rng);

        let rho1 = rho((x.0, x.1));
        let rho2 = rho((xprop.0, xprop.1));
        let ratio = rho2 / rho1;

        if rng.random_bool(ratio.min(1.0)) {
            count += 1;
            // Accept
            x = xprop;
        }

        if i % NSKIP == 0 {
            // Sample
            let mut fp = fp_mutex.lock().unwrap();
            fp.write(format!("{}\n", x.0 + SIZE * x.1).as_bytes()).unwrap();

            let r1 =
                ((x.0 as i32 - XMEAN as i32) * (x.0 as i32 - XMEAN as i32)
                 + (x.1 as i32 - YMEAN as i32) * (x.1 as i32 - YMEAN as i32)) as f64;
            let gauss = <f64>::exp(-r1/(2.0*SIGMA));
            let rho_x = rho((x.0, x.1));

            sample += gauss / rho_x;
        }
    }
    let _area_sq = NSAMP;
    let _area_circ = count;
    // A_s / A_c = pi r^2 / 4 r^2
    println!("count = {}", count);
    // A_g / A_c = pi / r^3
    let pi = (SIGMA / 2.0) * sample as f64 / ((NSAMP) as f64);
    println!("pi = {}", pi);
    pi
}

fn main() {
    let mut fp = File::create("distribution").unwrap();
    fp.write(format!("# {} {} {} ", SIZE, NSAMP, NWARM).as_bytes()).unwrap();
    let shared_fp = Arc::new(Mutex::new(fp));
    let mut threads = Vec::with_capacity(N_RUNS);
    let approx_pi = Arc::new(Mutex::new(0.0));
    for i in 0..N_RUNS {
        let fp = shared_fp.clone();
        let approx_pi_local = approx_pi.clone();
        let thread_handle = thread::spawn(move || {
            let pi = compute_pi(fp, i as u64);
            let mut handle = approx_pi_local.lock().unwrap();
            *handle += pi;
        });
        threads.push(thread_handle);
    }

    for thread in threads.into_iter() {
        thread.join().expect("Bruh");
    }
    let pi = approx_pi.lock().unwrap();
    println!("Approx pi = {}", *pi / N_RUNS as f64);
}
