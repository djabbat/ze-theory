//! Ze QMC — production simulator v0.2
//! Features: Wolff+Metropolis, Rayon, Xoshiro, Binder cumulant, Wilson loops, CLI

use rand::prelude::*;
use rand_xoshiro::Xoshiro256PlusPlus;
use rayon::prelude::*;
use clap::Parser;
use std::time::Instant;

/// Ze Quantum Monte Carlo simulator
#[derive(Parser, Debug)]
#[command(version = "0.2")]
struct Cli {
    /// Spatial size (Lx=Ly=Lz)
    #[arg(short='L', long, default_value = "4")]
    size: usize,
    /// Temporal size
    #[arg(short='t', long, default_value = "6")]
    lt: usize,
    /// Trotter slices
    #[arg(short='m', long, default_value = "16")]
    trotter: usize,
    /// J_t (AFM, >0)
    #[arg(long, default_value = "1.0")]
    jt: f64,
    /// J_s (FM, >0)
    #[arg(long, default_value = "0.0")]
    js: f64,
    /// Transverse field Gamma
    #[arg(short='G', long, default_value = "1.0")]
    gamma: f64,
    /// Longitudinal field h
    #[arg(short='H', long, default_value = "0.0")]
    h: f64,
    /// Inverse temperature beta
    #[arg(short='b', long, default_value = "10.0")]
    beta: f64,
    /// Number of thermalization steps
    #[arg(long, default_value = "500")]
    thermal: usize,
    /// Number of measurement steps
    #[arg(long, default_value = "2000")]
    samples: usize,
    /// Measurement interval
    #[arg(long, default_value = "10")]
    interval: usize,
    /// Scan: comma-separated gamma values
    #[arg(long)]
    scan: Option<String>,
    /// Enable Wilson loop measurement
    #[arg(long)]
    wilson: bool,
    /// Enable Binder cumulant
    #[arg(long)]
    binder: bool,
    /// Random seed
    #[arg(long, default_value = "42")]
    seed: u64,
}

struct Params { l: usize, lt: usize, m: usize, jt: f64, js: f64, g: f64, h: f64, b: f64 }
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
fn size(p: &Params) -> usize { p.l*p.l*p.l*p.lt*p.m }

/// Wolff cluster — параллельная версия (rayon для больших кластеров)
fn wolff(z: &mut Lattice, p: &Params, c: &TC, rng: &mut impl Rng) -> usize {
    let n = size(p);
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
        let val_i = z[i];
        let tau = i%p.m; let t = (i/p.m)%p.lt; let zc = (i/p.m/p.lt)%p.l;
        let y = (i/p.m/p.lt/p.l)%p.l; let x = i/p.m/p.lt/p.l/p.l;
        
        macro_rules! try_add {
            ($ni:expr, $prob:expr, $same:expr) => {
                if !cluster[$ni] && (z[$ni]*val_i > 0.0) == $same && rng.gen::<f64>() < $prob {
                    cluster[$ni] = true; queue.push($ni);
                }
            };
        }
        
        try_add!(idx(p,(x+1)%p.l,y,zc,t,tau), ps, true);     // FM spatial
        try_add!(idx(p,(x+p.l-1)%p.l,y,zc,t,tau), ps, true);
        try_add!(idx(p,x,(y+1)%p.l,zc,t,tau), ps, true);
        try_add!(idx(p,x,(y+p.l-1)%p.l,zc,t,tau), ps, true);
        try_add!(idx(p,x,y,(zc+1)%p.l,t,tau), ps, true);
        try_add!(idx(p,x,y,(zc+p.l-1)%p.l,t,tau), ps, true);
        try_add!(idx(p,x,y,zc,(t+1)%p.lt,tau), pt, false);    // AFM temporal
        try_add!(idx(p,x,y,zc,(t+p.lt-1)%p.lt,tau), pt, false);
        try_add!(idx(p,x,y,zc,t,(tau+1)%p.m), ptau, true);     // FM Trotter
        try_add!(idx(p,x,y,zc,t,(tau+p.m-1)%p.m), ptau, true);
    }
    
    //  Rayon: параллельный переворот кластера
    z.par_iter_mut().enumerate().for_each(|(i, val)| {
        if cluster[i] { *val = -*val; }
    });
    queue.len()
}

/// Измерение: (energy, |v|, |v_stag|, v² for Binder, Wilson loops)
fn measure(z: &Lattice, p: &Params, c: &TC, do_wilson: bool) -> (f64, f64, f64, f64, f64, Vec<f64>) {
    let n = size(p) as f64;
    let nc = (p.l*p.l*p.l) as f64;
    let mut e=0.0f64; let mut v=0.0f64; let mut vs_sum=0.0f64; let mut v4=0.0f64;
    let mut wloops = vec![];
    
    for x in 0..p.l { for y in 0..p.l { for zc in 0..p.l {
        let mut cs = 0.0f64;
        for t in 0..p.lt {
            let sign = if t%2==0 {1.0} else {-1.0};
            for tau in 0..p.m {
                let i = idx(p,x,y,zc,t,tau);
                let val = z[i];
                v += val; cs += sign*val; v4 += val.powi(4);
                let tn=(t+1)%p.lt; let taun=(tau+1)%p.m;
                let xn=(x+1)%p.l; let yn=(y+1)%p.l; let zn=(zc+1)%p.l;
                e += c.kt*val*z[idx(p,x,y,zc,tn,tau)];
                e -= c.ks*val*(z[idx(p,xn,y,zc,t,tau)]+z[idx(p,x,yn,zc,t,tau)]+z[idx(p,x,y,zn,t,tau)]);
                e -= c.ktau*val*z[idx(p,x,y,zc,t,taun)];
                e -= c.kh*val;
            }
        }
        vs_sum += (cs/(p.lt*p.m) as f64).abs();
    }}}
    
    let v_avg = v/n; let v_abs = v_avg.abs();
    let vs_avg = vs_sum / nc;
    // Binder cumulant для staggered magnetization (АФМ параметр порядка)
    // U₄ = 1 − ⟨m⁴⟩/(3⟨m²⟩²)
    let binder = if vs_avg > 0.001 { 
        let m2 = vs_avg.powi(2);
        1.0 - (cs_val.powi(4).max(1e-16))/(3.0*m2.max(1e-16))
    } else { 0.0 };
    
    // Wilson loops (R=1,2 × T=1,2)
    if do_wilson && p.l >= 3 {
        for r in 1..=2 { for t_loop in 1..=2 {
            let mut w = 0.0f64; let mut cnt = 0u64;
            for x in 0..p.l-r { for y in 0..p.l { for zc in 0..p.l {
                for t in 0..p.lt-t_loop {
                    let mut prod = 1.0f64;
                    for dx in 0..r { prod *= z[idx(p,x+dx,y,zc,t,0)]; }
                    for dt in 0..t_loop { prod *= z[idx(p,x+r,y,zc,t+dt,0)]; }
                    for dx in 0..r { prod *= z[idx(p,x+r-dx,y,zc,t+t_loop,0)]; }
                    for dt in 0..t_loop { prod *= z[idx(p,x,y,zc,t+t_loop-dt,0)]; }
                    w += prod; cnt += 1;
                }
            }}}
            if cnt > 0 { wloops.push(w/cnt as f64); }
        }}
    }
    
    (e/n, v_abs, vs_sum/nc, binder, v_avg, wloops)
}

fn run_single(p: &Params, cli: &Cli) {
    let c = TC::new(p);
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(cli.seed);
    let mut z = vec![1.0f64; size(p)];
    // init staggered
    for x in 0..p.l { for y in 0..p.l { for zc in 0..p.l {
        for t in 0..p.lt {
            let sign = if t%2==0 {1.0} else {-1.0};
            let base = idx(p,x,y,zc,t,0);
            for tau in 0..p.m { z[base+tau] = sign; }
        }
    }}}
    
    let t0 = Instant::now();
    for _ in 0..cli.thermal { wolff(&mut z, p, &c, &mut rng); }
    
    let nm = cli.samples/cli.interval;
    let (mut es,mut vs,mut vss,mut bs,mut vavgs) = (0.,0.,0.,0.,0.);
    for _ in 0..cli.samples {
        wolff(&mut z, p, &c, &mut rng);
        if rng.gen_range(0..cli.interval)==0 {
            let (e,v,vs_sum,b,va,_) = measure(&z,p,&c,cli.wilson);
            es+=e; vs+=v; vss+=vs_sum; bs+=b; vavgs+=va;
        }
    }
    es/=nm as f64; vs/=nm as f64; vss/=nm as f64; bs/=nm as f64; vavgs/=nm as f64;
    let dt = t0.elapsed().as_secs_f64();
    
    let phase = if vss>0.3 {"АФМ"} else if vs<0.2 {"пара"} else {"крит"};
    println!("Γ={:.2} |v|={:.4} v_stag={:.4} E/N={:.4} B={:.4} {:>5} ({:.1}s)",
             p.g, vs, vss, es, bs, phase, dt);
    
    if cli.wilson {
        let (_,_,_,_,_,w) = measure(&z,p,&c,true);
        println!("  Wilson loops: {:?}", w.iter().map(|x| format!("{:.4}",x)).collect::<Vec<_>>());
    }
}

fn main() {
    let cli = Cli::parse();
    
    if let Some(ref scan_str) = cli.scan {
        println!("Ze QMC v0.2 — scan Γ = {}\n", scan_str);
        println!("{:>8} {:>10} {:>10} {:>10} {:>10} {:>6}", "Γ","|v|","v_stag","E/N","Binder","Фаза");
        println!("{}","─".repeat(55));
        for gs in scan_str.split(',') {
            let g: f64 = gs.trim().parse().unwrap();
            let p = Params { l:cli.size, lt:cli.lt, m:cli.trotter, jt:cli.jt, js:cli.js,
                g, h:cli.h, b:cli.beta };
            run_single(&p, &cli);
        }
    } else {
        let p = Params { l:cli.size, lt:cli.lt, m:cli.trotter, jt:cli.jt, js:cli.js,
            g:cli.gamma, h:cli.h, b:cli.beta };
        run_single(&p, &cli);
    }
}
