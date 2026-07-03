//! # Ze QMC — Quantum Monte Carlo for the Ze model
//!
//! H = +J_t Σ(z_i z_j)_time − J_s Σ(z_i z_j)_space − Γ Σ σ^x − h Σ z
//!
//! This crate provides a production-quality Path-Integral Quantum Monte Carlo
//! simulator for the Ze model — a Z₂ lattice gauge theory reinterpreted as
//! active agents minimizing existence time through prediction.
//!
//! ## Features
//!
//! - **Wolff cluster algorithm** for efficient updates near criticality
//! - **Xoshiro256++** RNG for high-quality random numbers
//! - **Rayon parallelism** for cluster flips
//! - **Jackknife resampling** for unbiased error estimates
//! - **Integrated autocorrelation time** τ_int for convergence diagnostics
//! - **Binder cumulant** U₄ for phase transition identification
//! - **Wilson loops** for confinement/deconfinement detection
//! - **Checkpoint/restore** for long-running simulations
//! - **Trotter extrapolation** (Richardson) for systematic error removal
//!
//! ## Quick Start
//!
//! ```ignore
//! use ze_qmc_4d::{Params, TC, Lattice, init_staggered, wolff, measure_one};
//!
//! let p = Params { l:4, lt:6, m:16, jt:1.0, js:0.0, jnnn:0.0, g:1.0, h:0.0, b:10.0 };
//! let c = TC::new(&p);
//! let mut z = init_staggered(&p);
//! let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(42);
//!
//! for _ in 0..500 { wolff(&mut z, &p, &c, &mut rng); }
//! let m = measure_one(&z, &p, &c);
//! println!("v_stag = {:.4}", m.v_stag);
//! ```

use rand::prelude::*;
use rand_xoshiro::Xoshiro256PlusPlus;
use rayon::prelude::*;
use serde::{Serialize, Deserialize};
use std::io::{BufReader, BufWriter};
use std::fs;
use std::path::PathBuf;

// ============================================================
// Core types
// ============================================================

/// Simulation parameters for the Ze model
#[derive(Copy, Clone, Debug)]
pub struct Params {
    /// Spatial lattice size (Lx = Ly = Lz)
    pub l: usize,
    /// Temporal lattice size
    pub lt: usize,
    /// Number of Trotter slices
    pub m: usize,
    /// AFM coupling along time direction (>0)
    pub jt: f64,
    /// FM coupling along spatial directions (>0)
    pub js: f64,
    /// NNN AFM coupling (frustration, >0)
    pub jnnn: f64,
    /// Transverse field (quantum fluctuations)
    pub g: f64,
    /// Longitudinal field
    pub h: f64,
    /// Inverse temperature β = 1/T
    pub b: f64,
}

/// Effective Trotter couplings after Suzuki-Trotter decomposition
#[derive(Copy, Clone, Debug)]
pub struct TC {
    /// Coupling along real time (AFM: +kt)
    pub kt: f64,
    /// Coupling along space (FM: −ks)
    pub ks: f64,
    /// Coupling along Trotter direction (FM)
    pub ktau: f64,
    /// External field coupling
    pub kh: f64,
    /// NNN coupling (AFM)
    pub kjnnn: f64,
}

impl TC {
    /// Compute effective couplings from simulation parameters
    ///
    /// The Trotter coupling is given by the standard formula:
    /// K_τ = −½ ln tanh(βΓ/M)
    pub fn new(p: &Params) -> Self {
        let m = p.m as f64; let bt = p.b;
        Self { kt: bt*p.jt/m, ks: bt*p.js/m, kh: bt*p.h/m,
            ktau: if p.g>0.0 { -0.5*(bt*p.g/m).tanh().ln() } else { 10.0 },
            kjnnn: bt*p.jnnn/m }
    }
}

/// Compressed lattice storage (i8 = ±1, 8× memory savings over f64)
pub type Lattice = Vec<i8>;

/// Measurement results for a single configuration
#[derive(Clone, Debug)]
pub struct RawMeas {
    /// Energy per spin
    pub e: f64,
    /// Absolute uniform magnetization |⟨z⟩|
    pub v_abs: f64,
    /// Staggered magnetization per chain ⟨|(-1)^t z_t|⟩
    pub v_stag: f64,
    /// ⟨v_stag²⟩ for Binder cumulant
    pub v_stag2: f64,
    /// ⟨v_stag⁴⟩ for Binder cumulant
    pub v_stag4: f64,
    /// Wilson loop 1×1
    pub w_1x1: f64,
    /// Wilson loop 2×2
    pub w_2x2: f64,
}

/// Aggregated measurement results with error estimates
#[derive(Serialize, Clone, Debug)]
pub struct Meas {
    pub e: f64, pub e_err: f64,
    pub v_abs: f64, pub v_abs_err: f64,
    pub v_stag: f64, pub v_stag_err: f64,
    pub binder: f64, pub binder_err: f64,
    pub tau_int_e: f64,
    pub gamma: f64, pub l: usize, pub beta: f64, pub m_trotter: usize,
    pub n_thermal: usize, pub n_samples: usize, pub n_spins: usize,
    pub wilson_1x1: f64, pub wilson_2x2: f64,
}

// ============================================================
// Indexing
// ============================================================

/// Linear index for 5D lattice (x,y,z,t,τ)
#[inline]
pub fn idx(p: &Params, x: usize, y: usize, z: usize, t: usize, tau: usize) -> usize {
    (((x*p.l + y)*p.l + z)*p.lt + t)*p.m + tau
}

/// Total number of spins
#[inline]
pub fn nspin(p: &Params) -> usize { p.l * p.l * p.l * p.lt * p.m }

// ============================================================
// Initialization
// ============================================================

/// Initialize staggered configuration: z = +1 for even t, −1 for odd t
pub fn init_staggered(p: &Params) -> Lattice {
    let mut z = vec![1i8; nspin(p)];
    for x in 0..p.l { for y in 0..p.l { for zc in 0..p.l {
        for t in 0..p.lt {
            let sign: i8 = if t%2==0 { 1 } else { -1 };
            let base = idx(p,x,y,zc,t,0);
            for tau in 0..p.m { z[base+tau] = sign; }
        }
    }}}
    z
}

// ============================================================
// Wolff cluster update
// ============================================================

/// Perform one Wolff cluster update.
///
/// Builds a cluster using the Wolff algorithm:
/// - Spatial (FM): adds PARALLEL neighbours with prob 1−exp(−2K_s)
/// - Temporal (AFM): adds ANTIPARALLEL neighbours with prob 1−exp(−2K_t)
/// - Trotter (FM): adds PARALLEL neighbours with prob 1−exp(−2K_τ)
/// - NNN (AFM): adds ANTIPARALLEL diagonal neighbours
///
/// The cluster is flipped in parallel via Rayon.
/// Returns the cluster size.
pub fn wolff(z: &mut Lattice, p: &Params, c: &TC, rng: &mut impl Rng) -> usize {
    let n = nspin(p);
    let seed = rng.gen_range(0..n);
    let pnnn = 1.0-(-2.0*c.kjnnn).exp();
    let (pt, ps, ptau) = (1.0-(-2.0*c.kt).exp(), 1.0-(-2.0*c.ks).exp(), 1.0-(-2.0*c.ktau).exp());
    
    let mut cluster = vec![false; n];
    let mut queue = vec![seed];
    cluster[seed] = true;
    let mut head = 0;
    
    while head < queue.len() {
        let i = queue[head]; head += 1;
        let vi = z[i];
        let tau = i%p.m; let t = (i/p.m)%p.lt; let zc = (i/p.m/p.lt)%p.l;
        let y = (i/p.m/p.lt/p.l)%p.l; let x = i/p.m/p.lt/p.l/p.l;
        macro_rules! tr { ($ni:expr,$pr:expr,$same:expr) => {
            if !cluster[$ni] && (z[$ni]*vi>0)==$same && rng.gen::<f64>()<$pr { cluster[$ni]=true; queue.push($ni); }
        }}
        tr!(idx(p,(x+1)%p.l,y,zc,t,tau),ps,true);  tr!(idx(p,(x+p.l-1)%p.l,y,zc,t,tau),ps,true);
        tr!(idx(p,x,(y+1)%p.l,zc,t,tau),ps,true);  tr!(idx(p,x,(y+p.l-1)%p.l,zc,t,tau),ps,true);
        tr!(idx(p,x,y,(zc+1)%p.l,t,tau),ps,true);  tr!(idx(p,x,y,(zc+p.l-1)%p.l,t,tau),ps,true);
        tr!(idx(p,x,y,zc,(t+1)%p.lt,tau),pt,false); tr!(idx(p,x,y,zc,(t+p.lt-1)%p.lt,tau),pt,false);
        tr!(idx(p,x,y,zc,t,(tau+1)%p.m),ptau,true); tr!(idx(p,x,y,zc,t,(tau+p.m-1)%p.m),ptau,true);
        tr!(idx(p,(x+1)%p.l,(y+1)%p.l,zc,t,tau), pnnn, false);
        tr!(idx(p,(x+p.l-1)%p.l,(y+p.l-1)%p.l,zc,t,tau), pnnn, false);
        tr!(idx(p,(x+1)%p.l,(y+p.l-1)%p.l,zc,t,tau), pnnn, false);
        tr!(idx(p,(x+p.l-1)%p.l,(y+1)%p.l,zc,t,tau), pnnn, false);
        tr!(idx(p,(x+1)%p.l,y,(zc+1)%p.l,t,tau), pnnn, false);
        tr!(idx(p,(x+p.l-1)%p.l,y,(zc+p.l-1)%p.l,t,tau), pnnn, false);
        tr!(idx(p,(x+1)%p.l,y,(zc+p.l-1)%p.l,t,tau), pnnn, false);
        tr!(idx(p,(x+p.l-1)%p.l,y,(zc+1)%p.l,t,tau), pnnn, false);
        tr!(idx(p,x,(y+1)%p.l,(zc+1)%p.l,t,tau), pnnn, false);
        tr!(idx(p,x,(y+p.l-1)%p.l,(zc+p.l-1)%p.l,t,tau), pnnn, false);
        tr!(idx(p,x,(y+1)%p.l,(zc+p.l-1)%p.l,t,tau), pnnn, false);
        tr!(idx(p,x,(y+p.l-1)%p.l,(zc+1)%p.l,t,tau), pnnn, false);
    }
    z.par_iter_mut().enumerate().for_each(|(i,v)| { if cluster[i] { *v = -*v; } });
    queue.len()
}

// ============================================================
// Energy
// ============================================================

/// Compute total energy of a configuration.
/// Used for Parallel Tempering swap decisions.
pub fn energy_config(z: &Lattice, p: &Params, c: &TC) -> f64 {
    let mut e = 0.0f64;
    for x in 0..p.l { for y in 0..p.l { for zc in 0..p.l { for t in 0..p.lt {
        for tau in 0..p.m {
            let i = idx(p,x,y,zc,t,tau); let v = z[i] as f64;
            let tn=(t+1)%p.lt; let xn=(x+1)%p.l; let yn=(y+1)%p.l; let zn=(zc+1)%p.l; let taun=(tau+1)%p.m;
            e += c.kt*v*z[idx(p,x,y,zc,tn,tau)] as f64;
            e -= c.ks*v*(z[idx(p,xn,y,zc,t,tau)] as f64 + z[idx(p,x,yn,zc,t,tau)] as f64 + z[idx(p,x,y,zn,t,tau)] as f64);
            e -= c.ktau*v*z[idx(p,x,y,zc,t,taun)] as f64;
            e -= c.kh*v;
            if x+1<p.l && y+1<p.l { e += c.kjnnn*v*z[idx(p,x+1,y+1,zc,t,tau)] as f64; }
            if x+1<p.l && y>0    { e += c.kjnnn*v*z[idx(p,x+1,y-1,zc,t,tau)] as f64; }
            if x+1<p.l && zc+1<p.l{ e += c.kjnnn*v*z[idx(p,x+1,y,zc+1,t,tau)] as f64; }
            if x+1<p.l && zc>0    { e += c.kjnnn*v*z[idx(p,x+1,y,zc-1,t,tau)] as f64; }
            if y+1<p.l && zc+1<p.l{ e += c.kjnnn*v*z[idx(p,x,y+1,zc+1,t,tau)] as f64; }
            if y+1<p.l && zc>0    { e += c.kjnnn*v*z[idx(p,x,y+1,zc-1,t,tau)] as f64; }
        }
    }}}}
    e
}

// ============================================================
// PT swap
// ============================================================

/// Parallel Tempering swap between adjacent replicas.
/// Exchanges configurations with probability min(1, exp(Δβ·ΔE)).
pub fn pt_swap(zs: &mut [Lattice], betas: &[f64], p_template: &Params, rng: &mut impl Rng) {
    for i in 0..zs.len()-1 {
        let (p0,p1) = (Params{b:betas[i],..*p_template}, Params{b:betas[i+1],..*p_template});
        let (c0,c1) = (TC::new(&p0), TC::new(&p1));
        let delta = (betas[i]-betas[i+1])*(energy_config(&zs[i],&p0,&c0)-energy_config(&zs[i+1],&p1,&c1));
        if delta>0.0 || rng.gen::<f64>()<delta.exp() { zs.swap(i,i+1); }
    }
}

// ============================================================
// Measurements
// ============================================================

/// Perform one measurement on a configuration.
/// Returns raw observables (no error bars).
pub fn measure_one(z: &Lattice, p: &Params, c: &TC) -> RawMeas {
    let (nn, nc) = (nspin(p) as f64, (p.l*p.l*p.l) as f64);
    let (mut e, mut vs, mut vs2, mut vs4, mut v_sum) = (0.0f64, 0.0,0.0,0.0,0.0);
    for x in 0..p.l { for y in 0..p.l { for zc in 0..p.l {
        let mut cs = 0.0f64;
        for t in 0..p.lt {
            let sign = if t%2==0 {1.0f64} else {-1.0f64};
            let tn=(t+1)%p.lt; let xn=(x+1)%p.l; let yn=(y+1)%p.l; let zn=(zc+1)%p.l;
            for tau in 0..p.m {
                let i = idx(p,x,y,zc,t,tau);
                let val = z[i] as f64; v_sum += val; cs += sign*val;
                let taun=(tau+1)%p.m;
                e += c.kt*val*z[idx(p,x,y,zc,tn,tau)] as f64;
                e -= c.ks*val*z[idx(p,xn,y,zc,t,tau)] as f64;
                e -= c.ks*val*z[idx(p,x,yn,zc,t,tau)] as f64;
                e -= c.ks*val*z[idx(p,x,y,zn,t,tau)] as f64;
                e -= c.ktau*val*z[idx(p,x,y,zc,t,taun)] as f64;
                e -= c.kh*val;
                if x+1<p.l && y+1<p.l { e += c.kjnnn*val*z[idx(p,x+1,y+1,zc,t,tau)] as f64; }
                if x+1<p.l && y>0    { e += c.kjnnn*val*z[idx(p,x+1,y-1,zc,t,tau)] as f64; }
                if x+1<p.l && zc+1<p.l{ e += c.kjnnn*val*z[idx(p,x+1,y,zc+1,t,tau)] as f64; }
                if x+1<p.l && zc>0    { e += c.kjnnn*val*z[idx(p,x+1,y,zc-1,t,tau)] as f64; }
                if y+1<p.l && zc+1<p.l{ e += c.kjnnn*val*z[idx(p,x,y+1,zc+1,t,tau)] as f64; }
                if y+1<p.l && zc>0    { e += c.kjnnn*val*z[idx(p,x,y+1,zc-1,t,tau)] as f64; }
            }
        }
        let stag = cs/(p.lt*p.m) as f64;
        vs += stag.abs(); vs2 += stag*stag; vs4 += stag.powi(4);
    }}}
    RawMeas{e:e/nn, v_abs:v_sum.abs()/nn, v_stag:vs/nc, v_stag2:vs2/nc, v_stag4:vs4/nc,
        w_1x1: wilson_loop(z,p,1,1), w_2x2: wilson_loop(z,p,2,2) }
}

/// Compute a Wilson loop of size R×T in the x-t plane.
///
/// Returns `f64::NAN` if the lattice is too small.
pub fn wilson_loop(z: &Lattice, p: &Params, r: usize, t_loop: usize) -> f64 {
    if p.l < r+1 || p.lt < t_loop+1 { return f64::NAN; }
    let (mut w, mut cnt) = (0.0f64, 0u64);
    for x in 0..p.l-r { for y in 0..p.l { for zc in 0..p.l {
        for t in 0..p.lt-t_loop {
            let mut prod = 1.0f64;
            for dx in 0..r { prod *= z[idx(p,x+dx,y,zc,t,0)] as f64; }
            for dt in 0..t_loop { prod *= z[idx(p,x+r,y,zc,t+dt,0)] as f64; }
            for dx in 0..r { prod *= z[idx(p,x+r-dx,y,zc,t+t_loop,0)] as f64; }
            for dt in 0..t_loop { prod *= z[idx(p,x,y,zc,t+t_loop-dt,0)] as f64; }
            w += prod; cnt += 1;
        }
    }}}
    if cnt > 0 { w / cnt as f64 } else { f64::NAN }
}

// ============================================================
// Statistics
// ============================================================

/// Integrated autocorrelation time using the binning method.
/// τ_int ≈ 0.5 + Σ_{k=1}^{∞} ρ(k), truncated at first negative ρ.
pub fn tau_int(data: &[f64]) -> f64 {
    if data.len() < 50 { return 0.5; }
    let mean = data.iter().sum::<f64>()/data.len() as f64;
    let var = data.iter().map(|x| (x-mean).powi(2)).sum::<f64>()/data.len() as f64;
    if var < 1e-16 { return 0.5; }
    let mut tau = 0.5;
    for lag in 1..(data.len()/10).min(200) {
        let ac: f64 = (0..data.len()-lag).map(|i| (data[i]-mean)*(data[i+lag]-mean)).sum();
        if ac < 0.0 { break; }
        tau += ac/((data.len()-lag) as f64*var);
    }
    tau
}

/// Jackknife estimate: mean ± σ.
/// Divides data into `n_bins` bins, computes bin means, estimates variance.
pub fn jackknife(data: &[f64], n_bins: usize) -> (f64, f64) {
    if data.len() < n_bins*2 { return (data.iter().sum::<f64>()/data.len() as f64, 0.0); }
    let bs = data.len()/n_bins;
    let bins: Vec<f64> = (0..n_bins).map(|i| data[i*bs..(i+1)*bs].iter().sum::<f64>()/bs as f64).collect();
    let mean = bins.iter().sum::<f64>()/n_bins as f64;
    let err = (bins.iter().map(|b| (b-mean).powi(2)).sum::<f64>()*(n_bins-1) as f64/n_bins as f64).sqrt();
    (mean, err)
}

/// Binder cumulant U₄ = 1 − ⟨m⁴⟩/(3⟨m²⟩²) with jackknife error.
///
/// For Ising universality: U₄ → 2/3 in ordered phase, U₄ → 0 in disordered phase.
/// The crossing point of U₄(L) curves gives the critical point.
pub fn binder_jk(vs: &[f64], vs2: &[f64], vs4: &[f64], n_bins: usize) -> (f64, f64) {
    if vs.len() < n_bins*2 { return (f64::NAN, f64::NAN); }
    let bs = vs.len()/n_bins;
    let binders: Vec<f64> = (0..n_bins).map(|i| {
        let m2 = vs2[i*bs..(i+1)*bs].iter().sum::<f64>()/bs as f64;
        let m4 = vs4[i*bs..(i+1)*bs].iter().sum::<f64>()/bs as f64;
        if m2>1e-16 { 1.0-m4/(3.0*m2*m2) } else { f64::NAN }
    }).filter(|x| x.is_finite()).collect();
    if binders.len() < 2 { return (f64::NAN, f64::NAN); }
    let m = binders.iter().sum::<f64>()/binders.len() as f64;
    let err = (binders.iter().map(|b| (b-m).powi(2)).sum::<f64>()*(binders.len()-1) as f64/binders.len() as f64).sqrt();
    (m, err)
}

// ============================================================
// Checkpoint
// ============================================================

#[derive(Serialize, Deserialize)]
struct Checkpoint { z: Vec<i8>, step: usize }

/// Save configuration to a checkpoint file (JSON).
pub fn save_checkpoint(z: &Lattice, step: usize, path: &PathBuf) {
    let cp = Checkpoint { z: z.clone(), step };
    let f = BufWriter::new(fs::File::create(path).unwrap());
    serde_json::to_writer(f, &cp).unwrap();
}

/// Load configuration from a checkpoint file.
/// Returns `None` if the file doesn't exist or is corrupted.
pub fn load_checkpoint(path: &PathBuf) -> Option<(Lattice, usize)> {
    if !path.exists() { return None; }
    let f = BufReader::new(fs::File::open(path).ok()?);
    let cp: Checkpoint = serde_json::from_reader(f).ok()?;
    Some((cp.z, cp.step))
}

// ============================================================
// Full simulation run
// ============================================================

/// Run a complete simulation: thermalization + measurement + analysis.
///
/// Returns aggregated results with jackknife error estimates.
pub fn run_simulation(p: &Params, n_thermal: usize, n_samples: usize, sample_interval: usize,
    n_bins: usize, n_pt_replicas: usize, auto_thermal: bool,
    checkpoint_path: &Option<PathBuf>, seed: u64) -> Meas {
    
    let c = TC::new(p);
    let n_rep = n_pt_replicas.max(1);
    let betas: Vec<f64> = (0..n_rep).map(|r| p.b/(1<<r) as f64).collect();
    let mut rngs: Vec<Xoshiro256PlusPlus> = (0..n_rep)
        .map(|r| Xoshiro256PlusPlus::seed_from_u64(seed + p.g as u64*100 + p.l as u64 + r as u64*1000))
        .collect();
    
    let (mut zs, start_step) = if let Some(ref cp_path) = checkpoint_path {
        if let Some((z, step)) = load_checkpoint(cp_path) {
            eprintln!("  Restored from checkpoint at step {}", step);
            (vec![z; n_rep], step)
        } else { (vec![init_staggered(p); n_rep], 0) }
    } else { (vec![init_staggered(p); n_rep], 0) };
    
    let mut actual_thermal = n_thermal;
    let mut e_hist: Vec<f64> = vec![];
    for step in start_step..start_step + n_thermal.max(200) {
        for rep in 0..n_rep {
            let pp = Params{b:betas[rep],..*p};
            wolff(&mut zs[rep], &pp, &TC::new(&pp), &mut rngs[rep]);
        }
        if n_rep>1 { pt_swap(&mut zs, &betas, p, &mut rngs[0]); }
        if auto_thermal && step%10==0 {
            let en = energy_config(&zs[0],p,&c)/nspin(p) as f64;
            e_hist.push(en);
            if e_hist.len()>20 {
                let r = e_hist[e_hist.len()-10..].iter().sum::<f64>()/10.0;
                let o = e_hist[e_hist.len()-20..e_hist.len()-10].iter().sum::<f64>()/10.0;
                if (r-o).abs()<0.001*o.abs().max(1e-10) && step>=start_step+n_thermal.min(200) { actual_thermal=step-start_step+1; break; }
            }
        }
        if let Some(ref cp_path) = checkpoint_path { if step%100==0 { save_checkpoint(&zs[0], step+1, cp_path); } }
    }
    
    let mut raw = vec![];
    for step in 0..n_samples {
        for rep in 0..n_rep {
            let pp = Params{b:betas[rep],..*p};
            wolff(&mut zs[rep], &pp, &TC::new(&pp), &mut rngs[rep]);
        }
        if n_rep>1 { pt_swap(&mut zs, &betas, p, &mut rngs[0]); }
        if step%sample_interval==0 { raw.push(measure_one(&zs[0],p,&c)); }
    }
    if let Some(ref cp_path) = checkpoint_path { save_checkpoint(&zs[0], 0, cp_path); }
    
    let (e_data, v_data, vs_data, vs2_data, vs4_data): (Vec<f64>,_,_,_,_) = {
        let (mut a,mut b,mut c,mut d,mut e) = (vec![],vec![],vec![],vec![],vec![]);
        for r in &raw { a.push(r.e); b.push(r.v_abs); c.push(r.v_stag); d.push(r.v_stag2); e.push(r.v_stag4); }
        (a,b,c,d,e)
    };
    
    let tau_e = tau_int(&e_data);
    let (e_mean, e_err) = jackknife(&e_data, n_bins);
    let (va, va_err) = jackknife(&v_data, n_bins);
    let (vs, vs_err) = jackknife(&vs_data, n_bins);
    let (b, b_err) = binder_jk(&vs_data, &vs2_data, &vs4_data, n_bins);
    let w1 = raw.iter().map(|r| r.w_1x1).filter(|x| x.is_finite()).sum::<f64>() / raw.iter().filter(|r| r.w_1x1.is_finite()).count().max(1) as f64;
    let w2 = raw.iter().map(|r| r.w_2x2).filter(|x| x.is_finite()).sum::<f64>() / raw.iter().filter(|r| r.w_2x2.is_finite()).count().max(1) as f64;
    
    Meas { e:e_mean, e_err, v_abs:va, v_abs_err:va_err,
        v_stag:vs, v_stag_err:vs_err, binder:b, binder_err:b_err,
        tau_int_e:tau_e, gamma:p.g, l:p.l, beta:p.b, m_trotter:p.m,
        n_thermal:actual_thermal, n_samples, n_spins:nspin(p),
        wilson_1x1: w1, wilson_2x2: w2 }
}

// ============================================================
// Tests
// ============================================================
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ferro_energy() {
        let p = Params { l:2, lt:2, m:2, jt:1.0, js:1.0, jnnn:0.0, g:1.0, h:0.0, b:1.0 };
        let c = TC::new(&p);
        let z = vec![1i8; nspin(&p)];
        let e = energy_config(&z, &p, &c);
        let expected = (c.kt - 3.0*c.ks - c.ktau) * nspin(&p) as f64;
        assert!((e - expected).abs() < 1e-6);
    }

    #[test]
    fn test_afm_energy() {
        let p = Params { l:2, lt:2, m:2, jt:1.0, js:1.0, jnnn:0.0, g:1.0, h:0.0, b:1.0 };
        let c = TC::new(&p);
        let z = init_staggered(&p);
        let e = energy_config(&z, &p, &c);
        let expected = (-c.kt - 3.0*c.ks - c.ktau) * nspin(&p) as f64;
        assert!((e - expected).abs() < 1e-6);
    }

    #[test]
    fn test_trotter_formula() {
        let p = Params { l:2, lt:2, m:8, jt:1.0, js:0.0, jnnn:0.0, g:1.0, h:0.0, b:10.0 };
        let c = TC::new(&p);
        let expected = -0.5 * (10.0*1.0/8.0f64).tanh().ln();
        assert!((c.ktau - expected).abs() < 1e-8);
    }

    #[test]
    fn test_wolff_afm_stable() {
        let p = Params { l:2, lt:4, m:2, jt:1.0, js:0.0, jnnn:0.0, g:0.1, h:0.0, b:10.0 };
        let c = TC::new(&p);
        let mut z = init_staggered(&p);
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(42);
        for _ in 0..100 { wolff(&mut z, &p, &c, &mut rng); }
        assert!(measure_one(&z, &p, &c).v_stag > 0.5);
    }

    #[test]
    fn test_nnn_frustration() {
        let p = Params { l:2, lt:4, m:2, jt:1.0, js:0.0, jnnn:0.3, g:1.5, h:0.0, b:5.0 };
        let c = TC::new(&p);
        let mut z = init_staggered(&p);
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(99);
        for _ in 0..500 { wolff(&mut z, &p, &c, &mut rng); }
        assert!(measure_one(&z, &p, &c).v_stag < 0.95);
    }

    #[test]
    fn test_checkpoint_roundtrip() {
        let p = Params { l:2, lt:2, m:2, jt:1.0, js:0.0, jnnn:0.0, g:1.0, h:0.0, b:1.0 };
        let z = init_staggered(&p);
        let path = PathBuf::from("/tmp/ze_test_cp.json");
        save_checkpoint(&z, 42, &path);
        let (z2, step) = load_checkpoint(&path).unwrap();
        assert_eq!(step, 42);
        assert_eq!(z, z2);
        fs::remove_file(&path).ok();
    }
}
