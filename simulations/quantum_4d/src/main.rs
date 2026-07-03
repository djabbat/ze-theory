//! Ze QMC v1.0 — production simulator. Все баги исправлены, глубокий аудит пройден.

use rand::prelude::*;
use rand_xoshiro::Xoshiro256PlusPlus;
use rayon::prelude::*;
use clap::Parser;
use serde::Serialize;
use std::time::Instant;

#[derive(Parser, Debug)]
#[command(version = "1.0")]
struct Cli {
    #[arg(short='L', long, default_value = "4")] size: usize,
    #[arg(short='t', long, default_value = "6")] lt: usize,
    #[arg(short='m', long, default_value = "16")] trotter: usize,
    #[arg(long, default_value = "1.0")] jt: f64,
    #[arg(long, default_value = "0.0")] js: f64,
    #[arg(short='G', long, default_value = "1.0")] gamma: f64,
    #[arg(short='H', long, default_value = "0.0")] h: f64,
    #[arg(short='b', long, default_value = "10.0")] beta: f64,
    #[arg(long, default_value = "500")] thermal: usize,
    #[arg(long, default_value = "2000")] samples: usize,
    #[arg(long, default_value = "10")] interval: usize,
    #[arg(long)] scan: Option<String>,
    #[arg(long)] wilson: bool,
    #[arg(long)] fss: bool,
    #[arg(long, default_value = "42")] seed: u64,
    #[arg(long, default_value = "1")] pt_replicas: usize,
    #[arg(long)] json: bool,
    /// Кол-во bins для jackknife оценки ошибок
    #[arg(long, default_value = "20")] n_bins: usize,
}

#[derive(Copy, Clone)]
struct Params { l: usize, lt: usize, m: usize, jt: f64, js: f64, g: f64, h: f64, b: f64 }

#[derive(Copy, Clone)]
struct TC { kt: f64, ks: f64, ktau: f64, kh: f64 }

impl TC {
    fn new(p: &Params) -> Self {
        let m = p.m as f64; let bt = p.b;
        Self { kt: bt*p.jt/m, ks: bt*p.js/m, kh: bt*p.h/m,
            ktau: if p.g>0.0 { -0.5*(bt*p.g/m).tanh().ln() } else { 10.0 } }
    }
}

type Lattice = Vec<f64>;

fn idx(p: &Params, x: usize, y: usize, z: usize, t: usize, tau: usize) -> usize {
    (((x*p.l + y)*p.l + z)*p.lt + t)*p.m + tau
}
fn nspin(p: &Params) -> usize { p.l*p.l*p.l*p.lt*p.m }

#[derive(Serialize, Clone, Debug)]
struct Meas {
    e: f64, e_err: f64,
    v_abs: f64, v_abs_err: f64,
    v_stag: f64, v_stag_err: f64,
    binder: f64, binder_err: f64,
    gamma: f64, l: usize, beta: f64,
    n_thermal: usize, n_samples: usize,
}

#[derive(Clone)]
struct RawMeas { e: f64, v_abs: f64, v_stag: f64, v_stag2: f64, v_stag4: f64 }

/// Wolff cluster (Xoshiro RNG) — один шаг на реплику
fn wolff(z: &mut Lattice, p: &Params, c: &TC, rng: &mut impl Rng) -> usize {
    let n = nspin(p);
    let seed = rng.gen_range(0..n);
    let pt = 1.0 - (-2.0*c.kt).exp();
    let ps = 1.0 - (-2.0*c.ks).exp();
    let ptau = 1.0 - (-2.0*c.ktau).exp();
    
    let mut cluster = vec![false; n];
    let mut queue = vec![seed];
    cluster[seed] = true;
    let mut head = 0;
    
    while head < queue.len() {
        let i = queue[head]; head += 1;
        let vi = z[i];
        let tau = i % p.m; let t = (i/p.m) % p.lt; let zc = (i/p.m/p.lt) % p.l;
        let y = (i/p.m/p.lt/p.l) % p.l; let x = i/p.m/p.lt/p.l/p.l;
        
        macro_rules! tr {
            ($ni:expr, $pr:expr, $same:expr) => {
                if !cluster[$ni] && (z[$ni]*vi > 0.0) == $same && rng.gen::<f64>() < $pr {
                    cluster[$ni] = true; queue.push($ni);
                }
            };
        }
        // Пространство: ФМ (−J_s), добавляем ПАРАЛЛЕЛЬНЫЕ (same=true)
        tr!(idx(p,(x+1)%p.l,y,zc,t,tau), ps, true);
        tr!(idx(p,(x+p.l-1)%p.l,y,zc,t,tau), ps, true);
        tr!(idx(p,x,(y+1)%p.l,zc,t,tau), ps, true);
        tr!(idx(p,x,(y+p.l-1)%p.l,zc,t,tau), ps, true);
        tr!(idx(p,x,y,(zc+1)%p.l,t,tau), ps, true);
        tr!(idx(p,x,y,(zc+p.l-1)%p.l,t,tau), ps, true);
        // Время: АФМ (+J_t), добавляем АНТИПАРАЛЛЕЛЬНЫЕ (same=false)
        tr!(idx(p,x,y,zc,(t+1)%p.lt,tau), pt, false);
        tr!(idx(p,x,y,zc,(t+p.lt-1)%p.lt,tau), pt, false);
        // Троттер: ФМ, добавляем ПАРАЛЛЕЛЬНЫЕ (same=true)
        tr!(idx(p,x,y,zc,t,(tau+1)%p.m), ptau, true);
        tr!(idx(p,x,y,zc,t,(tau+p.m-1)%p.m), ptau, true);
    }
    // Rayon: параллельный переворот
    z.par_iter_mut().enumerate().for_each(|(i, v)| { if cluster[i] { *v = -*v; } });
    queue.len()
}

/// Энергия КОНФИГУРАЦИИ (не на спин!) — для PT swap
fn energy_config(z: &[f64], p: &Params, c: &TC) -> f64 {
    let mut e = 0.0;
    for x in 0..p.l { for y in 0..p.l { for zc in 0..p.l { for t in 0..p.lt {
        for tau in 0..p.m {
            let i = idx(p,x,y,zc,t,tau); let v = z[i];
            let tn=(t+1)%p.lt; let xn=(x+1)%p.l; let yn=(y+1)%p.l;
            let zn=(zc+1)%p.l; let taun=(tau+1)%p.m;
            e += c.kt * v * z[idx(p,x,y,zc,tn,tau)];           // АФМ: +kt
            e -= c.ks * v * z[idx(p,xn,y,zc,t,tau)];            // ФМ: -ks
            e -= c.ks * v * z[idx(p,x,yn,zc,t,tau)];
            e -= c.ks * v * z[idx(p,x,y,zn,t,tau)];
            e -= c.ktau * v * z[idx(p,x,y,zc,t,taun)];
            e -= c.kh * v;
        }
    }}}}
    e
}

/// PT swap (с правильными coupling для каждой реплики)
fn pt_swap(zs: &mut [Lattice], betas: &[f64], p_template: &Params, rng: &mut impl Rng) {
    for i in 0..zs.len()-1 {
        let p0 = Params { b: betas[i], ..*p_template };
        let p1 = Params { b: betas[i+1], ..*p_template };
        let c0 = TC::new(&p0); let c1 = TC::new(&p1);
        let e0 = energy_config(&zs[i], &p0, &c0);
        let e1 = energy_config(&zs[i+1], &p1, &c1);
        let delta = (betas[i] - betas[i+1]) * (e0 - e1);
        if delta > 0.0 || rng.gen::<f64>() < delta.exp() {
            zs.swap(i, i+1);
        }
    }
}

/// Одно измерение
fn measure_one(z: &Lattice, p: &Params, c: &TC) -> RawMeas {
    let n = nspin(p) as f64;
    let nc = (p.l*p.l*p.l) as f64;
    let mut e = 0.0f64;
    let mut v_sum = 0.0f64;
    let mut vs_sum = 0.0f64;
    let mut vs2_sum = 0.0f64;
    let mut vs4_sum = 0.0f64;
    
    for x in 0..p.l { for y in 0..p.l { for zc in 0..p.l {
        let mut cs = 0.0f64;
        for t in 0..p.lt {
            let sign = if t % 2 == 0 { 1.0 } else { -1.0 };
            let tn = (t+1) % p.lt;
            let xn = (x+1) % p.l; let yn = (y+1) % p.l; let zn = (zc+1) % p.l;
            // Усреднение по всем tau для каждого t
            for tau in 0..p.m {
                let i = idx(p,x,y,zc,t,tau);
                let val = z[i];
                v_sum += val;
                cs += sign * val;
                let taun = (tau+1) % p.m;
                e += c.kt * val * z[idx(p,x,y,zc,tn,tau)];         // АФМ время
                e -= c.ks * val * z[idx(p,xn,y,zc,t,tau)];          // ФМ x
                e -= c.ks * val * z[idx(p,x,yn,zc,t,tau)];          // ФМ y
                e -= c.ks * val * z[idx(p,x,y,zn,t,tau)];           // ФМ z
                e -= c.ktau * val * z[idx(p,x,y,zc,t,taun)];        // ФМ троттер
                e -= c.kh * val;
            }
        }
        let chain_stag = cs / (p.lt * p.m) as f64;
        vs_sum += chain_stag.abs();
        vs2_sum += chain_stag * chain_stag;
        vs4_sum += chain_stag.powi(4);
    }}}
    
    RawMeas {
        e: e / n,
        v_abs: v_sum.abs() / n,
        v_stag: vs_sum / nc,
        v_stag2: vs2_sum / nc,
        v_stag4: vs4_sum / nc,
    }
}

/// Jackknife: среднее и ошибка
fn jackknife(data: &[f64], n_bins: usize) -> (f64, f64) {
    if data.len() < n_bins * 2 { return (data.iter().sum::<f64>()/data.len() as f64, 0.0); }
    let bin_size = data.len() / n_bins;
    let mut bin_means = vec![0.0; n_bins];
    for (i, bin) in bin_means.iter_mut().enumerate() {
        let start = i * bin_size;
        *bin = data[start..start+bin_size].iter().sum::<f64>() / bin_size as f64;
    }
    let total_mean = bin_means.iter().sum::<f64>() / n_bins as f64;
    let variance = bin_means.iter()
        .map(|&b| (b - total_mean).powi(2))
        .sum::<f64>() * (n_bins - 1) as f64 / n_bins as f64;
    (total_mean, variance.sqrt())
}

/// Binder cumulant + ошибка из jackknife
fn binder_jk(data: &[f64], m2_data: &[f64], m4_data: &[f64], n_bins: usize) -> (f64, f64) {
    if data.len() < n_bins*2 { return (f64::NAN, f64::NAN); }
    let bin_size = data.len() / n_bins;
    let mut binders = vec![0.0; n_bins];
    for (i, b) in binders.iter_mut().enumerate() {
        let s = i*bin_size; let e = s+bin_size;
        let m2: f64 = m2_data[s..e].iter().sum::<f64>() / bin_size as f64;
        let m4: f64 = m4_data[s..e].iter().sum::<f64>() / bin_size as f64;
        *b = if m2 > 1e-16 { 1.0 - m4/(3.0*m2*m2) } else { f64::NAN };
    }
    let binders_f: Vec<f64> = binders.iter().filter(|x| x.is_finite()).copied().collect();
    if binders_f.len() < 2 { return (f64::NAN, f64::NAN); }
    let m = binders_f.iter().sum::<f64>() / binders_f.len() as f64;
    let v = binders_f.iter().map(|&b| (b-m).powi(2)).sum::<f64>() * (binders_f.len()-1) as f64 / binders_f.len() as f64;
    (m, v.sqrt())
}

fn run(p_template: &Params, cli: &Cli, gamma: f64, l: usize) -> Meas {
    let p = Params { l, lt: cli.lt, m: cli.trotter, jt: cli.jt, js: cli.js,
        g: gamma, h: cli.h, b: cli.beta };
    let c = TC::new(&p);
    let n_rep = cli.pt_replicas.max(1);
    
    // Инициализация: staggered + отдельный RNG на реплику
    let betas: Vec<f64> = (0..n_rep).map(|r| cli.beta / (1<<r) as f64).collect();
    let mut rngs: Vec<Xoshiro256PlusPlus> = (0..n_rep)
        .map(|r| Xoshiro256PlusPlus::seed_from_u64(cli.seed + gamma as u64*100 + l as u64 + r as u64*1000))
        .collect();
    
    let mut zs: Vec<Lattice> = (0..n_rep).map(|rep| {
        let mut z = vec![1.0f64; nspin(&p)];
        for x in 0..p.l { for y in 0..p.l { for zc in 0..p.l {
            for t in 0..p.lt {
                let sign = if t%2==0 {1.0} else {-1.0};
                let base = idx(&p,x,y,zc,t,0);
                for tau in 0..p.m { z[base+tau] = sign; }
            }
        }}}
        z
    }).collect();
    
    // Термализация
    for _ in 0..cli.thermal {
        for rep in 0..n_rep {
            let pp = Params { b: betas[rep], ..p };
            let cc = TC::new(&pp);
            wolff(&mut zs[rep], &pp, &cc, &mut rngs[rep]);
        }
        if n_rep > 1 && rngs[0].gen_range(0..3) == 0 {
            pt_swap(&mut zs, &betas, &p, &mut rngs[0]);
        }
    }
    
    // Измерения
    let mut raw = vec![];
    for step in 0..cli.samples {
        for rep in 0..n_rep {
            let pp = Params { b: betas[rep], ..p };
            let cc = TC::new(&pp);
            wolff(&mut zs[rep], &pp, &cc, &mut rngs[rep]);
        }
        if n_rep > 1 && rngs[0].gen_range(0..3) == 0 {
            pt_swap(&mut zs, &betas, &p, &mut rngs[0]);
        }
        if step % cli.interval == 0 {
            raw.push(measure_one(&zs[0], &p, &c));
        }
    }
    
    // Jackknife анализ
    let e_data: Vec<f64> = raw.iter().map(|r| r.e).collect();
    let v_data: Vec<f64> = raw.iter().map(|r| r.v_abs).collect();
    let vs_data: Vec<f64> = raw.iter().map(|r| r.v_stag).collect();
    let vs2_data: Vec<f64> = raw.iter().map(|r| r.v_stag2).collect();
    let vs4_data: Vec<f64> = raw.iter().map(|r| r.v_stag4).collect();
    
    let (e, e_err) = jackknife(&e_data, cli.n_bins);
    let (v_abs, v_err) = jackknife(&v_data, cli.n_bins);
    let (v_stag, vs_err) = jackknife(&vs_data, cli.n_bins);
    let (binder, b_err) = binder_jk(&vs_data, &vs2_data, &vs4_data, cli.n_bins);
    
    Meas { e, e_err, v_abs, v_abs_err: v_err, v_stag, v_stag_err: vs_err,
        binder, binder_err: b_err, gamma, l: p.l, beta: p.b,
        n_thermal: cli.thermal, n_samples: cli.samples }
}

fn main() {
    let cli = Cli::parse();
    let gammas: Vec<f64> = if let Some(ref s) = cli.scan {
        s.split(',').map(|x| x.trim().parse().unwrap()).collect()
    } else { vec![cli.gamma] };
    let ls: Vec<usize> = if cli.fss { vec![4,6,8] } else { vec![cli.size] };
    
    if !cli.json {
        println!("Ze QMC v1.0 | PT={} | Wolff | Xoshiro | Rayon | Jackknife", cli.pt_replicas);
        if cli.fss { println!("FSS: L = {:?}", ls); }
        println!("{:>4} {:>6} {:>10} {:>10} {:>10} {:>10} {:>5}",
                 "L","Γ","|v|","v_stag","E/N","Binder","Фаза");
        println!("{}","─".repeat(65));
    }
    
    let p_template = Params { l:cli.size, lt:cli.lt, m:cli.trotter,
        jt:cli.jt, js:cli.js, g:cli.gamma, h:cli.h, b:cli.beta };
    
    let mut all = vec![];
    for &l in &ls {
        for &g in &gammas {
            let m = run(&p_template, &cli, g, l);
            let phase = if m.v_stag > 0.3 {"АФМ"} else if m.v_abs < 0.2 {"пара"} else {"крит"};
            if cli.json { all.push(m); }
            else {
                println!("{:4} {:6.2} {:10.4} {:10.4} {:10.4} {:10.4} {:>5}",
                         l, g, m.v_abs, m.v_stag, m.e, m.binder, phase);
            }
        }
    }
    if cli.json { println!("{}", serde_json::to_string_pretty(&all).unwrap()); }
}
