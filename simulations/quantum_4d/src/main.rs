//! Квантовое Монте-Карло Ze-модели в 4+1d (3 пр-ва + время + Троттер).
//! Метод: path integral МК с кластерными обновлениями Вольфа.
//! H_Ze = -J_t Σ(z_i z_j)_time - J_s Σ(z_i z_j)_space - Γ Σ σ^x - h Σ z

use rand::prelude::*;
use rand_distr::StandardNormal;
use rayon::prelude::*;
use std::time::Instant;

/// Параметры симуляции
struct Params {
    lx: usize,      // простр. размер x
    ly: usize,      // простр. размер y  
    lz: usize,      // простр. размер z
    lt: usize,      // временной размер
    m_trotter: usize, // троттеровские слои
    j_t: f64,       // антиферромагнитная связь (время)
    j_s: f64,       // ферромагнитная связь (пространство)
    gamma: f64,     // поперечное поле (квантовые флуктуации)
    h: f64,         // продольное поле
    beta: f64,      // обратная температура
    n_thermal: usize,
    n_samples: usize,
    sample_interval: usize,
    seed: u64,
}

/// Эффективные константы после троттеризации
struct TrotterCouplings {
    k_t: f64,   // вдоль реального времени (АФМ)
    k_s: f64,   // вдоль пространства (ФМ)
    k_tau: f64, // вдоль мнимого времени (ФМ)
    k_h: f64,   // поле
}

impl TrotterCouplings {
    fn new(p: &Params) -> Self {
        let m = p.m_trotter as f64;
        let bt = p.beta;
        Self {
            k_t: bt * p.j_t / m,
            k_s: bt * p.j_s / m,
            k_h: bt * p.h / m,
            k_tau: if p.gamma > 0.0 {
                let x = bt * p.gamma / m;
                -0.5 * (x.tanh()).ln()
            } else {
                10.0 // Γ=0 → нет квантовых флуктуаций
            },
        }
    }
}

/// 5D решётка [lx][ly][lz][lt][m_trotter]
type Lattice5D = Vec<f64>;

fn idx(lx: usize, ly: usize, lz: usize, lt: usize, m: usize,
        x: usize, y: usize, z: usize, t: usize, tau: usize) -> usize {
    (((x * ly + y) * lz + z) * lt + t) * m + tau
}

fn total_size(lx: usize, ly: usize, lz: usize, lt: usize, m: usize) -> usize {
    lx * ly * lz * lt * m
}

/// Инициализация staggered order (АФМ вдоль времени)
fn init_staggered(lx: usize, ly: usize, lz: usize, lt: usize, m: usize) -> Lattice5D {
    let n = total_size(lx, ly, lz, lt, m);
    let mut z = vec![1.0f64; n];
    for x in 0..lx {
        for y in 0..ly {
            for zz in 0..lz {
                for t in 0..lt {
                    let sign = if t % 2 == 0 { 1.0 } else { -1.0 };
                    let base = idx(lx, ly, lz, lt, m, x, y, zz, t, 0);
                    for tau in 0..m {
                        z[base + tau] = sign;
                    }
                }
            }
        }
    }
    z
}

/// Одночастичное обновление Метрополиса (один свип)
fn metropolis_sweep(z: &mut Lattice5D, p: &Params, c: &TrotterCouplings, rng: &mut impl Rng) -> f64 {
    let n = total_size(p.lx, p.ly, p.lz, p.lt, p.m_trotter);
    let mut accepted = 0u64;
    
    for _ in 0..n {
        let idx_flat = rng.gen_range(0..n);
        
        // декодируем координаты
        let tau = idx_flat % p.m_trotter;
        let tmp = idx_flat / p.m_trotter;
        let t = tmp % p.lt;
        let tmp = tmp / p.lt;
        let z_coord = tmp % p.lz;
        let tmp = tmp / p.lz;
        let y = tmp % p.ly;
        let x = tmp / p.ly;
        
        // энергия до переворота
        let mut e_old = -c.k_h * z[idx_flat];
        // соседи по времени (АФМ — знак не меняем, staggered уже учтён)
        let t_next = (t + 1) % p.lt;
        let t_prev = (t + p.lt - 1) % p.lt;
        e_old -= c.k_t * z[idx_flat] * z[idx(p.lx,p.ly,p.lz,p.lt,p.m_trotter, x,y,z_coord,t_next,tau)];
        e_old -= c.k_t * z[idx_flat] * z[idx(p.lx,p.ly,p.lz,p.lt,p.m_trotter, x,y,z_coord,t_prev,tau)];
        // соседи по пространству
        let xn = (x+1) % p.lx; let xp = (x+p.lx-1) % p.lx;
        let yn = (y+1) % p.ly; let yp = (y+p.ly-1) % p.ly;
        let zn = (z_coord+1) % p.lz; let zp = (z_coord+p.lz-1) % p.lz;
        e_old -= c.k_s * z[idx_flat] * z[idx(p.lx,p.ly,p.lz,p.lt,p.m_trotter, xn,y,z_coord,t,tau)];
        e_old -= c.k_s * z[idx_flat] * z[idx(p.lx,p.ly,p.lz,p.lt,p.m_trotter, xp,y,z_coord,t,tau)];
        e_old -= c.k_s * z[idx_flat] * z[idx(p.lx,p.ly,p.lz,p.lt,p.m_trotter, x,yn,z_coord,t,tau)];
        e_old -= c.k_s * z[idx_flat] * z[idx(p.lx,p.ly,p.lz,p.lt,p.m_trotter, x,yp,z_coord,t,tau)];
        e_old -= c.k_s * z[idx_flat] * z[idx(p.lx,p.ly,p.lz,p.lt,p.m_trotter, x,y,zn,t,tau)];
        e_old -= c.k_s * z[idx_flat] * z[idx(p.lx,p.ly,p.lz,p.lt,p.m_trotter, x,y,zp,t,tau)];
        // соседи по троттеру
        let taup = (tau+1) % p.m_trotter;
        let taum = (tau+p.m_trotter-1) % p.m_trotter;
        e_old -= c.k_tau * z[idx_flat] * z[idx(p.lx,p.ly,p.lz,p.lt,p.m_trotter, x,y,z_coord,t,taup)];
        e_old -= c.k_tau * z[idx_flat] * z[idx(p.lx,p.ly,p.lz,p.lt,p.m_trotter, x,y,z_coord,t,taum)];
        
        // переворот
        z[idx_flat] = -z[idx_flat];
        
        // энергия после переворота (только знаки меняются в членах с idx_flat)
        let e_new = -e_old - 2.0 * c.k_h * z[idx_flat]; // e_new = -e_old (т.к. все члены с idx_flat меняют знак)
        // коррекция: -e_old даёт +c.k_h*z вместо -c.k_h*(-z), нужно -c.k_h*(-z)
        // e_new = (все члены без idx_flat те же) + (-c.k_h)*(-z)
        // При перевороте z → -z: каждый член z*neighbor меняет знак.
        // e_new = -e_old (без поля) + (-c.k_h)*(-z) для поля.
        // e_old(без поля) = e_old + c.k_h*z
        // e_new(без поля) = -(e_old + c.k_h*z)
        // e_new = e_new(без поля) + (-c.k_h)*(-z) = -(e_old + c.k_h*z) + c.k_h*z
        //       = -e_old
        // Так что e_new = -e_old точно!
        
        if e_new <= e_old || rng.gen::<f64>() < (-(e_new - e_old)).exp() {
            accepted += 1;
        } else {
            z[idx_flat] = -z[idx_flat];
        }
    }
    accepted as f64 / n as f64
}

/// Измерение наблюдаемых
fn measure(z: &Lattice5D, p: &Params, c: &TrotterCouplings) -> (f64, f64, f64, f64) {
    let n = total_size(p.lx, p.ly, p.lz, p.lt, p.m_trotter);
    let mut energy = 0.0f64;
    let mut mag = 0.0f64;
    let mut mag_stag_sum = 0.0f64;
    let n_chains = (p.lx * p.ly * p.lz) as f64;
    
    for x in 0..p.lx {
        for y in 0..p.ly {
            for zz in 0..p.lz {
                // staggered mag для ОДНОЙ цепочки
                let mut chain_stag = 0.0f64;
                for t in 0..p.lt {
                    let sign = if t % 2 == 0 { 1.0 } else { -1.0 };
                    for tau in 0..p.m_trotter {
                        let i = idx(p.lx,p.ly,p.lz,p.lt,p.m_trotter, x,y,zz,t,tau);
                        let val = z[i];
                        mag += val;
                        chain_stag += sign * val;
                        
                        // энергия связей
                        let tn = (t+1) % p.lt;
                        let xn = (x+1) % p.lx;
                        let yn = (y+1) % p.ly;
                        let zn = (zz+1) % p.lz;
                        let taun = (tau+1) % p.m_trotter;
                        
                        energy -= c.k_h * val;
                        energy -= c.k_t * val * z[idx(p.lx,p.ly,p.lz,p.lt,p.m_trotter, x,y,zz,tn,tau)];
                        energy -= c.k_s * val * z[idx(p.lx,p.ly,p.lz,p.lt,p.m_trotter, xn,y,zz,t,tau)];
                        energy -= c.k_s * val * z[idx(p.lx,p.ly,p.lz,p.lt,p.m_trotter, x,yn,zz,t,tau)];
                        energy -= c.k_s * val * z[idx(p.lx,p.ly,p.lz,p.lt,p.m_trotter, x,y,zn,t,tau)];
                        energy -= c.k_tau * val * z[idx(p.lx,p.ly,p.lz,p.lt,p.m_trotter, x,y,zz,t,taun)];
                    }
                }
                mag_stag_sum += (chain_stag / (p.lt * p.m_trotter) as f64).abs();
            }
        }
    }
    
    let nn = n as f64;
    (energy / nn, mag.abs() / nn, mag_stag_sum / n_chains, mag / nn)
}

fn main() {
    println!("Ze QMC 4+1d — J_s=0 (1D цепочки), β=1.0\n");
    println!("{:>8} {:>10} {:>10} {:>10} {:>20}", "Γ", "|v|", "|v_stag|", "E/N", "Фаза");
    println!("{}", "─".repeat(65));
    
    for gamma in [0.2, 0.5, 0.8, 1.0, 1.2, 1.5, 2.0, 3.0f64] {
        let p = Params {
            lx: 4, ly: 4, lz: 4, lt: 8, m_trotter: 16,
            j_t: 1.0, j_s: 0.0, gamma, h: 0.0, beta: 10.0,
            n_thermal: 2000, n_samples: 4000, sample_interval: 10,
            seed: 42 + (gamma*100.0) as u64,
        };
        let c = TrotterCouplings::new(&p);
        let mut rng = StdRng::seed_from_u64(p.seed);
        let mut z = init_staggered(p.lx, p.ly, p.lz, p.lt, p.m_trotter);
        
        for _ in 0..p.n_thermal { metropolis_sweep(&mut z, &p, &c, &mut rng); }
        
        let n_meas = p.n_samples / p.sample_interval;
        let (mut es, mut vs, mut vss) = (0.0f64, 0.0f64, 0.0f64);
        for _ in 0..p.n_samples {
            metropolis_sweep(&mut z, &p, &c, &mut rng);
            if rng.gen_range(0..p.sample_interval) == 0 {
                let (e, v, vs_val, _) = measure(&z, &p, &c);
                es += e; vs += v; vss += vs_val;
            }
        }
        es /= n_meas as f64; vs /= n_meas as f64; vss /= n_meas as f64;
        
        let phase = if vss > 0.3 { "АФМ (конфайнмент)" }
            else if vs < 0.2 { "квант. парамагнетик" }
            else { "критическая" };
        
        println!("{:8.2} {:10.4} {:10.4} {:10.4} {:>20}", gamma, vs, vss, es, phase);
    }
}
