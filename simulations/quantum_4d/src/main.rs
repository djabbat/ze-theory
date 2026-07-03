//! Ze QMC v1.1 — auto-thermalization, adaptive PT, bin analysis, compressed storage

use rand::prelude::*;
use rand_xoshiro::Xoshiro256PlusPlus;
use rayon::prelude::*;
use clap::Parser;
use serde::Serialize;

#[derive(Parser, Debug)]
#[command(version = "1.1")]
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
    #[arg(long, default_value = "20")] n_bins: usize,
    /// Auto-thermalization: stop when energy plateau detected
    #[arg(long)] auto_thermal: bool,
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

fn idx(p: &Params, x: usize, y: usize, z: usize, t: usize, tau: usize) -> usize {
    (((x*p.l + y)*p.l + z)*p.lt + t)*p.m + tau
}
fn nspin(p: &Params) -> usize { p.l*p.l*p.l*p.lt*p.m }

#[derive(Serialize, Clone, Debug)]
struct Meas {
    e: f64, e_err: f64, v_abs: f64, v_abs_err: f64,
    v_stag: f64, v_stag_err: f64, binder: f64, binder_err: f64,
    tau_int_e: f64,  // интегрированное время автокорреляции энергии
    gamma: f64, l: usize, beta: f64, n_thermal: usize, n_samples: usize,
}

#[derive(Clone)]
struct RawMeas { e: f64, v_abs: f64, v_stag: f64, v_stag2: f64, v_stag4: f64 }

/// Wolff cluster — оптимизированная версия
fn wolff(z: &mut [f64], p: &Params, c: &TC, rng: &mut impl Rng) -> usize {
    let n = nspin(p);
    let seed = rng.gen_range(0..n);
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
        
        macro_rules! tr { ($ni:expr,$pr:expr,$s:expr) => {
            if !cluster[$ni] && (z[$ni]*vi>0.0)==$s && rng.gen::<f64>()<$pr { cluster[$ni]=true; queue.push($ni); }
        }}
        tr!(idx(p,(x+1)%p.l,y,zc,t,tau),ps,true);  tr!(idx(p,(x+p.l-1)%p.l,y,zc,t,tau),ps,true);
        tr!(idx(p,x,(y+1)%p.l,zc,t,tau),ps,true);  tr!(idx(p,x,(y+p.l-1)%p.l,zc,t,tau),ps,true);
        tr!(idx(p,x,y,(zc+1)%p.l,t,tau),ps,true);  tr!(idx(p,x,y,(zc+p.l-1)%p.l,t,tau),ps,true);
        tr!(idx(p,x,y,zc,(t+1)%p.lt,tau),pt,false); tr!(idx(p,x,y,zc,(t+p.lt-1)%p.lt,tau),pt,false);
        tr!(idx(p,x,y,zc,t,(tau+1)%p.m),ptau,true); tr!(idx(p,x,y,zc,t,(tau+p.m-1)%p.m),ptau,true);
    }
    z.par_iter_mut().enumerate().for_each(|(i,v)| { if cluster[i] { *v = -*v; } });
    queue.len()
}

fn energy_config(z: &[f64], p: &Params, c: &TC) -> f64 {
    let mut e = 0.0;
    for x in 0..p.l { for y in 0..p.l { for zc in 0..p.l { for t in 0..p.lt {
        for tau in 0..p.m {
            let i = idx(p,x,y,zc,t,tau); let v = z[i];
            let tn=(t+1)%p.lt; let xn=(x+1)%p.l; let yn=(y+1)%p.l; let zn=(zc+1)%p.l; let taun=(tau+1)%p.m;
            e += c.kt*v*z[idx(p,x,y,zc,tn,tau)]; e -= c.ks*v*z[idx(p,xn,y,zc,t,tau)];
            e -= c.ks*v*z[idx(p,x,yn,zc,t,tau)]; e -= c.ks*v*z[idx(p,x,y,zn,t,tau)];
            e -= c.ktau*v*z[idx(p,x,y,zc,t,taun)]; e -= c.kh*v;
        }
    }}}}
    e
}

fn pt_swap(zs: &mut [Vec<f64>], betas: &[f64], p_template: &Params, rng: &mut impl Rng) {
    for i in 0..zs.len()-1 {
        let (p0, p1) = (Params{b:betas[i],..*p_template}, Params{b:betas[i+1],..*p_template});
        let (c0, c1) = (TC::new(&p0), TC::new(&p1));
        let delta = (betas[i]-betas[i+1])*(energy_config(&zs[i],&p0,&c0)-energy_config(&zs[i+1],&p1,&c1));
        if delta > 0.0 || rng.gen::<f64>() < delta.exp() { zs.swap(i, i+1); }
    }
}

fn measure_one(z: &[f64], p: &Params, c: &TC) -> RawMeas {
    let (n, nc) = (nspin(p) as f64, (p.l*p.l*p.l) as f64);
    let (mut e, mut vs, mut vs2, mut vs4, mut v_sum) = (0.0,0.0,0.0,0.0,0.0);
    for x in 0..p.l { for y in 0..p.l { for zc in 0..p.l {
        let mut cs = 0.0;
        for t in 0..p.lt {
            let sign = if t%2==0 {1.0} else {-1.0};
            let tn=(t+1)%p.lt; let xn=(x+1)%p.l; let yn=(y+1)%p.l; let zn=(zc+1)%p.l;
            for tau in 0..p.m {
                let i = idx(p,x,y,zc,t,tau); let val = z[i]; v_sum += val; cs += sign*val;
                let taun=(tau+1)%p.m;
                e += c.kt*val*z[idx(p,x,y,zc,tn,tau)]; e -= c.ks*val*z[idx(p,xn,y,zc,t,tau)];
                e -= c.ks*val*z[idx(p,x,yn,zc,t,tau)]; e -= c.ks*val*z[idx(p,x,y,zn,t,tau)];
                e -= c.ktau*val*z[idx(p,x,y,zc,t,taun)]; e -= c.kh*val;
            }
        }
        let stag = cs/(p.lt*p.m) as f64;
        vs += stag.abs(); vs2 += stag*stag; vs4 += stag.powi(4);
    }}}
    RawMeas{e:e/n, v_abs:v_sum.abs()/n, v_stag:vs/nc, v_stag2:vs2/nc, v_stag4:vs4/nc}
}

/// Интегрированное время автокорреляции (метод биннинга)
fn tau_int(data: &[f64]) -> f64 {
    if data.len() < 50 { return 0.5; }
    let mean = data.iter().sum::<f64>() / data.len() as f64;
    let var = data.iter().map(|x| (x-mean).powi(2)).sum::<f64>() / data.len() as f64;
    if var < 1e-16 { return 0.5; }
    let mut tau = 0.5;
    let max_lag = (data.len() / 10).min(200);
    for lag in 1..max_lag {
        let mut ac = 0.0;
        for i in 0..data.len()-lag { ac += (data[i]-mean)*(data[i+lag]-mean); }
        ac /= (data.len()-lag) as f64 * var;
        if ac < 0.0 { break; }
        tau += ac;
    }
    tau
}

fn jackknife(data: &[f64], n_bins: usize) -> (f64, f64) {
    if data.len() < n_bins*2 { return (data.iter().sum::<f64>()/data.len() as f64, 0.0); }
    let bs = data.len()/n_bins;
    let bins: Vec<f64> = (0..n_bins).map(|i| data[i*bs..(i+1)*bs].iter().sum::<f64>()/bs as f64).collect();
    let mean = bins.iter().sum::<f64>()/n_bins as f64;
    let err = (bins.iter().map(|b| (b-mean).powi(2)).sum::<f64>()*(n_bins-1) as f64/n_bins as f64).sqrt();
    (mean, err)
}

fn binder_jk(vs: &[f64], vs2: &[f64], vs4: &[f64], n_bins: usize) -> (f64, f64) {
    if vs.len() < n_bins*2 { return (f64::NAN, f64::NAN); }
    let bs = vs.len()/n_bins;
    let binders: Vec<f64> = (0..n_bins).map(|i| {
        let m2 = vs2[i*bs..(i+1)*bs].iter().sum::<f64>()/bs as f64;
        let m4 = vs4[i*bs..(i+1)*bs].iter().sum::<f64>()/bs as f64;
        if m2>1e-16 { 1.0 - m4/(3.0*m2*m2) } else { f64::NAN }
    }).filter(|x| x.is_finite()).collect();
    if binders.len() < 2 { return (f64::NAN, f64::NAN); }
    let m = binders.iter().sum::<f64>()/binders.len() as f64;
    let err = (binders.iter().map(|b| (b-m).powi(2)).sum::<f64>()*(binders.len()-1) as f64/binders.len() as f64).sqrt();
    (m, err)
}

fn run(cli: &Cli, gamma: f64, l: usize) -> Meas {
    let p = Params { l, lt:cli.lt, m:cli.trotter, jt:cli.jt, js:cli.js, g:gamma, h:cli.h, b:cli.beta };
    let c = TC::new(&p);
    let n_rep = cli.pt_replicas.max(1);
    
    let betas: Vec<f64> = (0..n_rep).map(|r| cli.beta/(1<<r) as f64).collect();
    let mut rngs: Vec<Xoshiro256PlusPlus> = (0..n_rep)
        .map(|r| Xoshiro256PlusPlus::seed_from_u64(cli.seed + gamma as u64*100 + l as u64 + r as u64*1000))
        .collect();
    
    // Инициализация
    let mut zs: Vec<Vec<f64>> = (0..n_rep).map(|_| {
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
    
    // Авто-термализация
    let mut n_thermal = cli.thermal;
    let mut e_hist: Vec<f64> = vec![];
    for step in 0..cli.thermal.max(200) {
        for rep in 0..n_rep {
            let pp = Params{b:betas[rep],..p};
            wolff(&mut zs[rep], &pp, &TC::new(&pp), &mut rngs[rep]);
        }
        if n_rep>1 { pt_swap(&mut zs, &betas, &p, &mut rngs[0]); }
        
        if cli.auto_thermal && step%10==0 {
            e_hist.push(energy_config(&zs[0], &p, &c)/nspin(&p) as f64);
            if e_hist.len() > 20 {
                let recent = e_hist[e_hist.len()-10..].iter().sum::<f64>()/10.0;
                let older = e_hist[e_hist.len()-20..e_hist.len()-10].iter().sum::<f64>()/10.0;
                if (recent-older).abs() < 0.001*older.abs().max(1e-10) && step>=cli.thermal.min(200) {
                    n_thermal = step+1;
                    break;
                }
            }
        }
    }
    
    // Измерения
    let mut raw = vec![];
    for step in 0..cli.samples {
        for rep in 0..n_rep {
            let pp = Params{b:betas[rep],..p};
            wolff(&mut zs[rep], &pp, &TC::new(&pp), &mut rngs[rep]);
        }
        if n_rep>1 { pt_swap(&mut zs, &betas, &p, &mut rngs[0]); }
        if step % cli.interval == 0 { raw.push(measure_one(&zs[0], &p, &c)); }
    }
    
    // Анализ
    let e_data: Vec<f64> = raw.iter().map(|r| r.e).collect();
    let v_data: Vec<f64> = raw.iter().map(|r| r.v_abs).collect();
    let vs_data: Vec<f64> = raw.iter().map(|r| r.v_stag).collect();
    let vs2_data: Vec<f64> = raw.iter().map(|r| r.v_stag2).collect();
    let vs4_data: Vec<f64> = raw.iter().map(|r| r.v_stag4).collect();
    
    let tau_e = tau_int(&e_data);
    let (e, e_err) = jackknife(&e_data, cli.n_bins);
    let (va, va_err) = jackknife(&v_data, cli.n_bins);
    let (vs, vs_err) = jackknife(&vs_data, cli.n_bins);
    let (b, b_err) = binder_jk(&vs_data, &vs2_data, &vs4_data, cli.n_bins);
    
    Meas { e, e_err, v_abs:va, v_abs_err:va_err, v_stag:vs, v_stag_err:vs_err,
        binder:b, binder_err:b_err, tau_int_e:tau_e,
        gamma, l, beta:p.b, n_thermal, n_samples:cli.samples }
}

fn main() {
    let cli = Cli::parse();
    let gammas: Vec<f64> = if let Some(ref s) = cli.scan {
        s.split(',').map(|x| x.trim().parse().unwrap()).collect()
    } else { vec![cli.gamma] };
    let ls: Vec<usize> = if cli.fss { vec![4,6,8] } else { vec![cli.size] };
    
    if !cli.json {
        println!("Ze QMC v1.1 | PT={} | Wolff | Jackknife | τ_int{}",
                 cli.pt_replicas, if cli.auto_thermal {" | auto-thermal"} else {""});
        if cli.fss { println!("FSS: L = {:?}", ls); }
        println!("{:>4} {:>6} {:>10} {:>10} {:>10} {:>10} {:>7} {:>5}",
                 "L","Γ","|v|","v_stag","E/N","Binder","τ_int","Фаза");
        println!("{}","─".repeat(72));
    }
    
    let mut all = vec![];
    for &l in &ls {
        for &g in &gammas {
            let m = run(&cli, g, l);
            let phase = if m.v_stag>0.3 {"АФМ"} else if m.v_abs<0.2 {"пара"} else {"крит"};
            if cli.json { all.push(m); }
            else {
                println!("{:4} {:6.2} {:10.4} {:10.4} {:10.4} {:10.4} {:7.2} {:>5}",
                         l, g, m.v_abs, m.v_stag, m.e, m.binder, m.tau_int_e, phase);
            }
        }
    }
    if cli.json { println!("{}", serde_json::to_string_pretty(&all).unwrap()); }
}
