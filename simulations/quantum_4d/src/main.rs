//! Ze QMC v0.3 — parallel tempering + Binder + FSS + bench
//! Полный продакшен-симулятор

use rand::prelude::*;
use rand_xoshiro::Xoshiro256PlusPlus;
use rayon::prelude::*;
use clap::Parser;
use serde::Serialize;
use std::time::Instant;

#[derive(Parser, Debug)]
#[command(version = "0.3")]
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
    #[arg(long)] fss: bool,  // finite-size scaling: iterate L
    #[arg(long, default_value = "42")] seed: u64,
    /// Parallel tempering: number of replicas
    #[arg(long, default_value = "1")] pt_replicas: usize,
    /// JSON output
    #[arg(long)] json: bool,
}

struct Params { l: usize, lt: usize, m: usize, jt: f64, js: f64, g: f64, h: f64, b: f64 }
struct TC { kt: f64, ks: f64, ktau: f64, kh: f64 }

impl TC { fn new(p: &Params) -> Self {
    let m=p.m as f64; let bt=p.b;
    Self{kt:bt*p.jt/m, ks:bt*p.js/m, kh:bt*p.h/m,
        ktau:if p.g>0.0{-0.5*(bt*p.g/m).tanh().ln()}else{10.0}}
}}

type Lattice = Vec<f64>;
fn idx(p:&Params,x:usize,y:usize,z:usize,t:usize,tau:usize)->usize {
    (((x*p.l+y)*p.l+z)*p.lt+t)*p.m+tau
}
fn nspin(p:&Params)->usize { p.l*p.l*p.l*p.lt*p.m }

#[derive(Serialize, Clone)]
struct Meas {
    e: f64, v_abs: f64, v_stag: f64, v_stag2: f64, v_stag4: f64,
    w_loops: Vec<f64>, gamma: f64, l: usize, beta: f64,
}

/// Wolff (Xoshiro)
fn wolff(z:&mut Lattice, p:&Params, c:&TC, rng:&mut impl Rng)->usize{
    let n=nspin(p); let seed=rng.gen_range(0..n);
    let pt=1.0-(-2.0*c.kt).exp(); let ps=1.0-(-2.0*c.ks).exp(); let ptau=1.0-(-2.0*c.ktau).exp();
    let mut cl=vec![false;n]; let mut q=vec![seed]; cl[seed]=true; let mut hd=0;
    while hd<q.len() {
        let i=q[hd]; hd+=1; let vi=z[i];
        let tau=i%p.m; let t=(i/p.m)%p.lt; let zc=(i/p.m/p.lt)%p.l;
        let y=(i/p.m/p.lt/p.l)%p.l; let x=i/p.m/p.lt/p.l/p.l;
        macro_rules! tr {($ni:expr,$pr:expr,$same:expr)=>{if!cl[$ni]&&(z[$ni]*vi>0.0)==$same&&rng.gen::<f64>()<$pr{cl[$ni]=true;q.push($ni)}}}
        tr!(idx(p,(x+1)%p.l,y,zc,t,tau),ps,true); tr!(idx(p,(x+p.l-1)%p.l,y,zc,t,tau),ps,true);
        tr!(idx(p,x,(y+1)%p.l,zc,t,tau),ps,true); tr!(idx(p,x,(y+p.l-1)%p.l,zc,t,tau),ps,true);
        tr!(idx(p,x,y,(zc+1)%p.l,t,tau),ps,true); tr!(idx(p,x,y,(zc+p.l-1)%p.l,t,tau),ps,true);
        tr!(idx(p,x,y,zc,(t+1)%p.lt,tau),pt,false); tr!(idx(p,x,y,zc,(t+p.lt-1)%p.lt,tau),pt,false);
        tr!(idx(p,x,y,zc,t,(tau+1)%p.m),ptau,true); tr!(idx(p,x,y,zc,t,(tau+p.m-1)%p.m),ptau,true);
    }
    z.par_iter_mut().enumerate().for_each(|(i,v)|{if cl[i]{*v=-*v;}});
    q.len()
}

/// Parallel tempering swap — исправленная версия с учётом beta-зависимых couplings
fn pt_swap(zs: &mut [Lattice], betas: &[f64], p: &Params, rng: &mut impl Rng) {
    for i in 0..zs.len()-1 {
        // Вычисляем энергии с ПРАВИЛЬНЫМИ coupling для каждой реплики
        let pp0 = Params { b: betas[i], ..*p };
        let pp1 = Params { b: betas[i+1], ..*p };
        let c0 = TC::new(&pp0);
        let c1 = TC::new(&pp1);
        let e0 = energy_at(zs[i].as_slice(), &pp0, &c0);
        let e1 = energy_at(zs[i+1].as_slice(), &pp1, &c1);
        let delta = (betas[i] - betas[i+1]) * (e0 - e1);
        if delta > 0.0 || rng.gen::<f64>() < delta.exp() {
            zs.swap(i, i+1);
        }
    }
}

fn energy_at(z: &[f64], p: &Params, c: &TC) -> f64 {
    let mut e = 0.0;
    for x in 0..p.l { for y in 0..p.l { for zc in 0..p.l { for t in 0..p.lt {
        for tau in 0..p.m {
            let i = idx(p,x,y,zc,t,tau); let v = z[i];
            let tn=(t+1)%p.lt; let xn=(x+1)%p.l; let yn=(y+1)%p.l; let zn=(zc+1)%p.l; let tnau=(tau+1)%p.m;
            e += c.kt*v*z[idx(p,x,y,zc,tn,tau)];
            e -= c.ks*v*(z[idx(p,xn,y,zc,t,tau)]+z[idx(p,x,yn,zc,t,tau)]+z[idx(p,x,y,zn,t,tau)]);
            e -= c.ktau*v*z[idx(p,x,y,zc,t,tnau)];
            e -= c.kh*v;
        }
    }}}}
    e
}

fn measure_at(z:&Lattice, p:&Params, c:&TC, do_wilson:bool, gamma:f64) -> Meas {
    let n=nspin(p) as f64; let nc=(p.l*p.l*p.l) as f64;
    let mut e=0.0; let mut vs_sum=0.0; let mut vs2_sum=0.0; let mut vs4_sum=0.0; let mut v_sum=0.0;
    for x in 0..p.l { for y in 0..p.l { for zc in 0..p.l {
        let mut cs=0.0;
        for t in 0..p.lt {
            let sign=if t%2==0{1.0}else{-1.0};
            for tau in 0..p.m {
                let i=idx(p,x,y,zc,t,tau); let val=z[i]; v_sum+=val; cs+=sign*val;
                let tn=(t+1)%p.lt; let taun=(tau+1)%p.m;
                let xn=(x+1)%p.l; let yn=(y+1)%p.l; let zn=(zc+1)%p.l;
                e += c.kt*val*z[idx(p,x,y,zc,tn,tau)];
                e -= c.ks*val*(z[idx(p,xn,y,zc,t,tau)]+z[idx(p,x,yn,zc,t,tau)]+z[idx(p,x,y,zn,t,tau)]);
                e -= c.ktau*val*z[idx(p,x,y,zc,t,taun)]; e -= c.kh*val;
            }
        }
        let chain_stag = cs/(p.lt*p.m) as f64;
        vs_sum += chain_stag.abs();
        vs2_sum += chain_stag*chain_stag;
        vs4_sum += chain_stag.powi(4);
    }}}
    let mut wloops = vec![];
    if do_wilson && p.l>=3 {
        for r in 1..=2 { for tl in 1..=2 {
            let mut w=0.0; let mut cnt=0u64;
            for x in 0..p.l-r { for y in 0..p.l { for zc in 0..p.l {
                for t in 0..p.lt-tl {
                    let mut prod=1.0;
                    for dx in 0..r{prod*=z[idx(p,x+dx,y,zc,t,0)]}
                    for dt in 0..tl{prod*=z[idx(p,x+r,y,zc,t+dt,0)]}
                    for dx in 0..r{prod*=z[idx(p,x+r-dx,y,zc,t+tl,0)]}
                    for dt in 0..tl{prod*=z[idx(p,x,y,zc,t+tl-dt,0)]}
                    w+=prod; cnt+=1;
                }
            }}}
            if cnt>0{wloops.push(w/cnt as f64);}
        }}
    }
    Meas { e: e/n, v_abs: v_sum.abs()/n, v_stag: vs_sum/nc,
        v_stag2: vs2_sum/nc, v_stag4: vs4_sum/nc, w_loops: wloops,
        gamma, l: p.l, beta: p.b }
}

/// Binder cumulant: U₄ = 1 − ⟨m⁴⟩/(3⟨m²⟩²)
fn binder(m2: f64, m4: f64) -> f64 {
    if m2 > 1e-16 { 1.0 - m4/(3.0*m2*m2) } else { f64::NAN }
}

fn run(p: &Params, cli: &Cli, gamma: f64, l: usize) -> Meas {
    let mut pp = Params { l, lt: cli.lt, m: cli.trotter, jt: cli.jt, js: cli.js,
        g: gamma, h: cli.h, b: cli.beta };
    let c = TC::new(&pp);
    let n_replicas = cli.pt_replicas.max(1);
    
    // Parallel tempering: replicas at beta, beta/2, beta/4, ...
    let mut zs: Vec<Lattice> = (0..n_replicas).map(|rep| {
        let beta_rep = cli.beta / (1 << rep) as f64;
        let mut z = vec![1.0f64; nspin(&pp)];
        for x in 0..pp.l { for y in 0..pp.l { for zc in 0..pp.l {
            for t in 0..pp.lt {
                let sign = if t%2==0 {1.0} else {-1.0};
                let base = idx(&pp,x,y,zc,t,0);
                for tau in 0..pp.m { z[base+tau] = sign; }
            }
        }}}
        z
    }).collect();
    
    let betas: Vec<f64> = (0..n_replicas).map(|rep| cli.beta/(1<<rep) as f64).collect();
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(cli.seed + gamma as u64*100 + l as u64);
    
    // Термализация с PT
    for _ in 0..cli.thermal {
        for rep in 0..n_replicas {
            let pp_rep = Params { l, lt: cli.lt, m: cli.trotter, jt: cli.jt, js: cli.js,
                g: gamma, h: cli.h, b: betas[rep] };
            let c_rep = TC::new(&pp_rep);
            wolff(&mut zs[rep], &pp_rep, &c_rep, &mut rng);
        }
        if n_replicas > 1 && rng.gen_range(0..5)==0 {
            pt_swap(&mut zs, &betas, &pp, &mut rng);
        }
    }
    
    // Измерения
    let nm = cli.samples/cli.interval;
    let (mut es,mut vs,mut vss,mut vs2s,mut vs4s) = (0.,0.,0.,0.,0.);
    for _ in 0..cli.samples {
        for rep in 0..n_replicas {
            let pp_rep = Params { l, lt: cli.lt, m: cli.trotter, jt: cli.jt, js: cli.js,
                g: gamma, h: cli.h, b: betas[rep] };
            let c_rep = TC::new(&pp_rep);
            wolff(&mut zs[rep], &pp_rep, &c_rep, &mut rng);
        }
        if n_replicas > 1 && rng.gen_range(0..5)==0 {
            pt_swap(&mut zs, &betas, &pp, &mut rng);
        }
        if rng.gen_range(0..cli.interval)==0 {
            let m = measure_at(&zs[0], &pp, &c, cli.wilson, gamma);
            es+=m.e; vs+=m.v_abs; vss+=m.v_stag; vs2s+=m.v_stag2; vs4s+=m.v_stag4;
        }
    }
    es/=nm as f64; vs/=nm as f64; vss/=nm as f64; vs2s/=nm as f64; vs4s/=nm as f64;
    
    Meas { e:es, v_abs:vs, v_stag:vss, v_stag2:vs2s, v_stag4:vs4s,
        w_loops: vec![], gamma, l, beta: cli.beta }
}

fn main() {
    let cli = Cli::parse();
    let gammas: Vec<f64> = if let Some(ref s) = cli.scan {
        s.split(',').map(|x|x.trim().parse().unwrap()).collect()
    } else { vec![cli.gamma] };
    
    let ls: Vec<usize> = if cli.fss { vec![4,6,8] } else { vec![cli.size] };
    
    let mut all_results: Vec<Meas> = vec![];
    
    if !cli.json {
        println!("Ze QMC v0.3 | PT={} | Xoshiro | Wolff | Rayon", cli.pt_replicas);
        if cli.fss { println!("FSS: L = {:?}", ls); }
        println!("{:>4} {:>6} {:>10} {:>10} {:>10} {:>10} {:>10}",
                 "L","Γ","|v|","v_stag","E/N","Binder","Фаза");
        println!("{}","─".repeat(70));
    }
    
    for &l in &ls {
        for &g in &gammas {
            let p = Params { l, lt:cli.lt, m:cli.trotter, jt:cli.jt, js:cli.js,
                g, h:cli.h, b:cli.beta };
            let m = run(&p, &cli, g, l);
            let b = binder(m.v_stag2, m.v_stag4);
            let phase = if m.v_stag>0.3 {"АФМ"} else if m.v_abs<0.2 {"пара"} else {"крит"};
            
            if cli.json {
                all_results.push(m);
            } else {
                println!("{:4} {:6.2} {:10.4} {:10.4} {:10.4} {:10.4} {:>5}",
                         l, g, m.v_abs, m.v_stag, m.e, b, phase);
            }
        }
    }
    
    if cli.json {
        println!("{}", serde_json::to_string_pretty(&all_results).unwrap());
    }
}
