//! Ze QMC 4+1d — аудит-версия с Wolff кластерами и правильным знаком J_t
//! H = +J_t Σ(z_i z_j)_time − J_s Σ(z_i z_j)_space − Γ Σ σ^x − h Σ z

use rand::prelude::*;
use std::time::Instant;

struct Params { lx: usize, ly: usize, lz: usize, lt: usize, m: usize,
    j_t: f64, j_s: f64, gamma: f64, h: f64, beta: f64,
    n_thermal: usize, n_samples: usize, sample_interval: usize, seed: u64 }

struct TC { k_t: f64, k_s: f64, k_tau: f64, k_h: f64 }

impl TC {
    fn new(p: &Params) -> Self {
        let m = p.m as f64; let bt = p.beta;
        Self { k_t: bt * p.j_t / m, k_s: bt * p.j_s / m, k_h: bt * p.h / m,
            k_tau: if p.gamma > 0.0 { -0.5 * (bt * p.gamma / m).tanh().ln() } else { 10.0 } }
    }
}

type Lattice = Vec<f64>;

fn idx(p: &Params, x: usize, y: usize, z: usize, t: usize, tau: usize) -> usize {
    (((x * p.ly + y) * p.lz + z) * p.lt + t) * p.m + tau
}
fn total(p: &Params) -> usize { p.lx * p.ly * p.lz * p.lt * p.m }

/// Кластер Wolff: строит кластер из seed, переворачивает его.
/// Возвращает размер кластера.
fn wolff_flip(z: &mut Lattice, p: &Params, c: &TC, rng: &mut impl Rng) -> usize {
    let n = total(p);
    let seed = rng.gen_range(0..n);
    
    // вероятности добавления
    let p_t = 1.0 - (-2.0 * c.k_t).exp();   // АФМ: +J_t → K_t>0, проверяем АНТИпараллельность
    let p_s = 1.0 - (-2.0 * c.k_s).exp();   // ФМ: −J_s → wait, k_s=βJ_s/M>0, связь −k_s*z_i*z_j
    // Для ФМ (−J_s): энергия = −k_s*z_i*z_j. Параллельные: −k_s. Антипараллельные: +k_s.
    // При перевороте одного спина: ΔE = 2*k_s*(z_i*z_j).
    // Кластер: добавляем если z_i == z_j (параллельные) с вероятностью 1−exp(−2*k_s)
    
    let mut cluster = vec![false; n];
    let mut queue = vec![seed];
    cluster[seed] = true;
    let mut head = 0;
    
    while head < queue.len() {
        let i = queue[head]; head += 1;
        let val_i = z[i];
        
        // Декодируем координаты
        let tau = i % p.m; let tmp = i / p.m;
        let t = tmp % p.lt; let tmp = tmp / p.lt;
        let zc = tmp % p.lz; let tmp = tmp / p.lz;
        let y = tmp % p.ly; let x = tmp / p.ly;
        
        // Соседи: (x±1,y,z,t,τ), (x,y±1,z,t,τ), (x,y,z±1,t,τ), (x,y,z,t±1,τ), (x,y,z,t,τ±1)
        let neighbors = [
            ((x+1)%p.lx, y, zc, t, tau, c.k_s, p_s, false),  // ФМ: + если равны
            ((x+p.lx-1)%p.lx, y, zc, t, tau, c.k_s, p_s, false),
            (x, (y+1)%p.ly, zc, t, tau, c.k_s, p_s, false),
            (x, (y+p.ly-1)%p.ly, zc, t, tau, c.k_s, p_s, false),
            (x, y, (zc+1)%p.lz, t, tau, c.k_s, p_s, false),
            (x, y, (zc+p.lz-1)%p.lz, t, tau, c.k_s, p_s, false),
            (x, y, zc, (t+1)%p.lt, tau, c.k_t, p_t, true),   // АФМ: + если НЕ равны
            (x, y, zc, (t+p.lt-1)%p.lt, tau, c.k_t, p_t, true),
            (x, y, zc, t, (tau+1)%p.m, c.k_tau, 1.0-( -2.0*c.k_tau).exp(), false),
            (x, y, zc, t, (tau+p.m-1)%p.m, c.k_tau, 1.0-( -2.0*c.k_tau).exp(), false),
        ];
        
        for (nx, ny, nz, nt, ntau, k_val, prob, is_afm) in &neighbors {
            let ni = idx(p, *nx, *ny, *nz, *nt, *ntau);
            if cluster[ni] { continue; }
            let val_j = z[ni];
            let same = (val_i * val_j) > 0.0;
            let should_add = if *is_afm { !same } else { same };
            if should_add && rng.gen::<f64>() < *prob {
                cluster[ni] = true;
                queue.push(ni);
            }
        }
    }
    
    // Переворот кластера
    for i in 0..n {
        if cluster[i] { z[i] = -z[i]; }
    }
    queue.len()
}

/// Измерение наблюдаемых
fn measure(z: &Lattice, p: &Params, c: &TC) -> (f64, f64, f64) {
    let n = total(p) as f64;
    let n_chains = (p.lx * p.ly * p.lz) as f64;
    let mut energy = 0.0f64;
    let mut v_sum = 0.0f64;
    let mut v_stag_sum = 0.0f64;
    
    for x in 0..p.lx { for y in 0..p.ly { for zc in 0..p.lz {
        let mut chain_stag = 0.0f64;
        for t in 0..p.lt {
            let sign = if t % 2 == 0 { 1.0 } else { -1.0 };
            for tau in 0..p.m {
                let i = idx(p, x, y, zc, t, tau);
                let val = z[i];
                v_sum += val;
                chain_stag += sign * val;
                let tn = (t+1) % p.lt; let taun = (tau+1) % p.m;
                let xn = (x+1) % p.lx; let yn = (y+1) % p.ly; let zn = (zc+1) % p.lz;
                energy += c.k_t * val * z[idx(p,x,y,zc,tn,tau)];      // АФМ: +k_t
                energy -= c.k_s * val * z[idx(p,xn,y,zc,t,tau)];       // ФМ: −k_s
                energy -= c.k_s * val * z[idx(p,x,yn,zc,t,tau)];
                energy -= c.k_s * val * z[idx(p,x,y,zn,t,tau)];
                energy -= c.k_tau * val * z[idx(p,x,y,zc,t,taun)];
                energy -= c.k_h * val;
            }
        }
        v_stag_sum += (chain_stag / (p.lt * p.m) as f64).abs();
    }}}
    (energy / n, v_sum.abs() / n, v_stag_sum / n_chains)
}

fn main() {
    println!("Ze QMC 4+1d АУДИТ — Wolff + правильный знак J_t\n");
    println!("{:>8} {:>10} {:>10} {:>10} {:>20}", "Γ", "|v|", "|v_stag|", "E/N", "Фаза");
    println!("{}", "─".repeat(65));
    
    for gamma in [0.2, 0.5, 0.8, 1.0, 1.5, 2.0, 3.0f64] {
        let p = Params { lx:4, ly:4, lz:4, lt:8, m:16,
            j_t:1.0, j_s:0.0, gamma, h:0.0, beta:10.0,
            n_thermal:500, n_samples:2000, sample_interval:10,
            seed:42 + (gamma*100.0) as u64 };
        let c = TC::new(&p);
        let mut rng = StdRng::seed_from_u64(p.seed);
        let mut z = vec![1.0f64; total(&p)];
        
        // Инициализация: staggered
        for x in 0..p.lx { for y in 0..p.ly { for zc in 0..p.lz {
            for t in 0..p.lt {
                let sign = if t % 2 == 0 { 1.0 } else { -1.0 };
                let base = idx(&p, x, y, zc, t, 0);
                for tau in 0..p.m { z[base + tau] = sign; }
            }
        }}}
        
        let t0 = Instant::now();
        for _ in 0..p.n_thermal { wolff_flip(&mut z, &p, &c, &mut rng); }
        
        let n_meas = p.n_samples / p.sample_interval;
        let (mut es, mut vs, mut vss) = (0.0, 0.0, 0.0);
        for step in 0..p.n_samples {
            wolff_flip(&mut z, &p, &c, &mut rng);
            if step % p.sample_interval == 0 {
                let (e, v_abs, v_stag) = measure(&z, &p, &c);
                es += e; vs += v_abs; vss += v_stag;
            }
        }
        es /= n_meas as f64; vs /= n_meas as f64; vss /= n_meas as f64;
        
        let phase = if vss > 0.3 { "АФМ (конфайнмент)" }
            else if vs < 0.2 { "квант. парамагнетик" } else { "критическая" };
        let dt = t0.elapsed().as_secs_f64();
        println!("{:8.2} {:10.4} {:10.4} {:10.4} {:>20}  ({:.1}s)", gamma, vs, vss, es, phase, dt);
    }
}
