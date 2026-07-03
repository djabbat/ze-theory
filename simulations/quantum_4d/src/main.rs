//! Ze QMC v1.2 — compressed i8 storage, Binder crossing, production-ready
//!
//! H = +J_t Σ(z_i z_j)_time − J_s Σ(z_i z_j)_space − Γ Σ σ^x − h Σ z
//!
//! Features: Wolff clusters, Xoshiro RNG, Rayon parallel, Jackknife errors,
//!           τ_int autocorrelation, auto-thermalization, PT swaps,
//!           Binder cumulant, FSS, Wilson loops, JSON output.
//!
//! Build: cargo build --release
//! Usage: ./ze-qmc-4d --scan "0.5,0.8,1.0,1.2,1.5" --fss --auto-thermal

use rand::prelude::*;
use rand_xoshiro::Xoshiro256PlusPlus;
use rayon::prelude::*;
use clap::Parser;
use serde::Serialize;

#[derive(Parser, Debug)]
#[command(version = "1.2")]
struct Cli {
    #[arg(short='L', long, default_value = "4")] size: usize,
    #[arg(short='t', long, default_value = "6")] lt: usize,
    #[arg(short='m', long, default_value = "16")] trotter: usize,
    #[arg(long, default_value = "1.0")] jt: f64,
    #[arg(long, default_value = "0.0")] js: f64,
    /// NNN AFM coupling (создаёт фрустрацию)
    #[arg(long, default_value = "0.0")] jnnn: f64,
    #[arg(short='G', long, default_value = "1.0")] gamma: f64,
    #[arg(short='H', long, default_value = "0.0")] h: f64,
    #[arg(short='b', long, default_value = "10.0")] beta: f64,
    #[arg(long, default_value = "500")] thermal: usize,
    #[arg(long, default_value = "2000")] samples: usize,
    #[arg(long, default_value = "10")] interval: usize,
    #[arg(long)] scan: Option<String>,
    #[arg(long)] fss: bool,
    #[arg(long, default_value = "42")] seed: u64,
    #[arg(long, default_value = "1")] pt_replicas: usize,
    #[arg(long)] json: bool,
    #[arg(long, default_value = "20")] n_bins: usize,
    #[arg(long)] auto_thermal: bool,
}

#[derive(Copy, Clone)]
struct Params { l: usize, lt: usize, m: usize, jt: f64, js: f64, jnnn: f64, g: f64, h: f64, b: f64 }

#[derive(Copy, Clone)]
struct TC { kt: f64, ks: f64, ktau: f64, kh: f64, kjnnn: f64 }

impl TC {
    fn new(p: &Params) -> Self {
        let m = p.m as f64; let bt = p.b;
        Self { kt: bt*p.jt/m, ks: bt*p.js/m, kh: bt*p.h/m,
            ktau: if p.g>0.0 { -0.5*(bt*p.g/m).tanh().ln() } else { 10.0 },
            kjnnn: bt*p.jnnn/m }
    }
}

// ============================================================
// Сжатое хранение: i8 вместо f64 (экономия 8×)
// ============================================================
/// Решётка в сжатом формате: каждый байт = ±1
/// Для Wolff и измерений распаковываем во временный буфер
type Lattice = Vec<i8>;

fn idx(p: &Params, x: usize, y: usize, z: usize, t: usize, tau: usize) -> usize {
    (((x*p.l + y)*p.l + z)*p.lt + t)*p.m + tau
}
fn nspin(p: &Params) -> usize { p.l * p.l * p.l * p.lt * p.m }

fn init_staggered(p: &Params) -> Lattice {
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
// Wolff cluster — работает напрямую с i8
// ============================================================
fn wolff(z: &mut Lattice, p: &Params, c: &TC, rng: &mut impl Rng) -> usize {
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
    // Rayon: параллельный flip
    z.par_iter_mut().enumerate().for_each(|(i,v)| { if cluster[i] { *v = -*v; } });
    queue.len()
}

// ============================================================
// Энергия конфигурации (для PT). Работает с i8, возвращает f64.
// ============================================================
fn energy_config(z: &Lattice, p: &Params, c: &TC) -> f64 {
    let mut e = 0.0f64;
    for x in 0..p.l { for y in 0..p.l { for zc in 0..p.l { for t in 0..p.lt {
        for tau in 0..p.m {
            let i = idx(p,x,y,zc,t,tau);
            let v = z[i] as f64;
            let tn=(t+1)%p.lt; let xn=(x+1)%p.l; let yn=(y+1)%p.l; let zn=(zc+1)%p.l; let taun=(tau+1)%p.m;
            e += c.kt*v*z[idx(p,x,y,zc,tn,tau)] as f64;
            e -= c.ks*v*(z[idx(p,xn,y,zc,t,tau)] as f64 + z[idx(p,x,yn,zc,t,tau)] as f64 + z[idx(p,x,y,zn,t,tau)] as f64);
            e -= c.ktau*v*z[idx(p,x,y,zc,t,taun)] as f64;
            e -= c.kh*v;
        }
    }}}}
    e
}

// ============================================================
// PT swap
// ============================================================
fn pt_swap(zs: &mut [Lattice], betas: &[f64], p_template: &Params, rng: &mut impl Rng) {
    for i in 0..zs.len()-1 {
        let (p0,p1) = (Params{b:betas[i],..*p_template}, Params{b:betas[i+1],..*p_template});
        let (c0,c1) = (TC::new(&p0), TC::new(&p1));
        let delta = (betas[i]-betas[i+1])*(energy_config(&zs[i],&p0,&c0)-energy_config(&zs[i+1],&p1,&c1));
        if delta>0.0 || rng.gen::<f64>()<delta.exp() { zs.swap(i,i+1); }
    }
}

// ============================================================
// Измерения
// ============================================================
#[derive(Clone)]
struct RawMeas { e: f64, v_abs: f64, v_stag: f64, v_stag2: f64, v_stag4: f64,
    w_1x1: f64, w_1x2: f64, w_2x1: f64, w_2x2: f64 }

fn measure_one(z: &Lattice, p: &Params, c: &TC) -> RawMeas {
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
            }
        }
        let stag = cs/(p.lt*p.m) as f64;
        vs += stag.abs(); vs2 += stag*stag; vs4 += stag.powi(4);
    }}}
    RawMeas{e:e/nn, v_abs:v_sum.abs()/nn, v_stag:vs/nc, v_stag2:vs2/nc, v_stag4:vs4/nc,
        w_1x1: wilson_loop(z,p,1,1), w_1x2: wilson_loop(z,p,1,2),
        w_2x1: wilson_loop(z,p,2,1), w_2x2: wilson_loop(z,p,2,2) }
}

/// Петля Вильсона R×T в x-t плоскости
fn wilson_loop(z: &Lattice, p: &Params, r: usize, t_loop: usize) -> f64 {
    if p.l < r+1 || p.lt < t_loop+1 { return f64::NAN; }
    let (mut w, mut cnt) = (0.0f64, 0u64);
    for x in 0..p.l-r {
        for y in 0..p.l { for zc in 0..p.l {
            for t in 0..p.lt-t_loop {
                let mut prod = 1.0f64;
                for dx in 0..r { prod *= z[idx(p,x+dx,y,zc,t,0)] as f64; }
                for dt in 0..t_loop { prod *= z[idx(p,x+r,y,zc,t+dt,0)] as f64; }
                for dx in 0..r { prod *= z[idx(p,x+r-dx,y,zc,t+t_loop,0)] as f64; }
                for dt in 0..t_loop { prod *= z[idx(p,x,y,zc,t+t_loop-dt,0)] as f64; }
                w += prod; cnt += 1;
            }
        }}
    }
    if cnt > 0 { w / cnt as f64 } else { f64::NAN }
}

// ============================================================
// Статистика
// ============================================================
fn tau_int(data: &[f64]) -> f64 {
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
        if m2>1e-16 { 1.0-m4/(3.0*m2*m2) } else { f64::NAN }
    }).filter(|x| x.is_finite()).collect();
    if binders.len() < 2 { return (f64::NAN, f64::NAN); }
    let m = binders.iter().sum::<f64>()/binders.len() as f64;
    let err = (binders.iter().map(|b| (b-m).powi(2)).sum::<f64>()*(binders.len()-1) as f64/binders.len() as f64).sqrt();
    (m, err)
}

// ============================================================
// Основной прогон
// ============================================================
#[derive(Serialize, Clone, Debug)]
struct Meas {
    e: f64, e_err: f64, v_abs: f64, v_abs_err: f64,
    v_stag: f64, v_stag_err: f64, binder: f64, binder_err: f64,
    tau_int_e: f64, gamma: f64, l: usize, beta: f64,
    n_thermal: usize, n_samples: usize, n_spins: usize,
    wilson_1x1: f64, wilson_1x2: f64, wilson_2x2: f64,
}

fn run(cli: &Cli, gamma: f64, l: usize) -> Meas {
    let p = Params { l, lt:cli.lt, m:cli.trotter, jt:cli.jt, js:cli.js, jnnn:cli.jnnn, g:gamma, h:cli.h, b:cli.beta };
    let c = TC::new(&p);
    let n_rep = cli.pt_replicas.max(1);
    let betas: Vec<f64> = (0..n_rep).map(|r| cli.beta/(1<<r) as f64).collect();
    let mut rngs: Vec<Xoshiro256PlusPlus> = (0..n_rep)
        .map(|r| Xoshiro256PlusPlus::seed_from_u64(cli.seed+gamma as u64*100+l as u64+r as u64*1000))
        .collect();
    
    let mut zs: Vec<Lattice> = (0..n_rep).map(|_| init_staggered(&p)).collect();
    
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
            let en = energy_config(&zs[0],&p,&c)/nspin(&p) as f64;
            e_hist.push(en);
            if e_hist.len()>20 {
                let r = e_hist[e_hist.len()-10..].iter().sum::<f64>()/10.0;
                let o = e_hist[e_hist.len()-20..e_hist.len()-10].iter().sum::<f64>()/10.0;
                if (r-o).abs()<0.001*o.abs().max(1e-10) && step>=cli.thermal.min(200) { n_thermal=step+1; break; }
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
        if step%cli.interval==0 { raw.push(measure_one(&zs[0],&p,&c)); }
    }
    
    let (e_data, v_data, vs_data, vs2_data, vs4_data): (Vec<f64>,_,_,_,_) = {
        let (mut a,mut b,mut c,mut d,mut e) = (vec![],vec![],vec![],vec![],vec![]);
        for r in &raw { a.push(r.e); b.push(r.v_abs); c.push(r.v_stag); d.push(r.v_stag2); e.push(r.v_stag4); }
        (a,b,c,d,e)
    };
    
    let tau_e = tau_int(&e_data);
    let (e_mean, e_err) = jackknife(&e_data, cli.n_bins);
    let (va, va_err) = jackknife(&v_data, cli.n_bins);
    let (vs, vs_err) = jackknife(&vs_data, cli.n_bins);
    let (b, b_err) = binder_jk(&vs_data, &vs2_data, &vs4_data, cli.n_bins);
    let w_data: Vec<f64> = raw.iter().map(|r| r.w_1x1).filter(|x| x.is_finite()).collect();
    let w_1x1 = if w_data.len()>0 { w_data.iter().sum::<f64>()/w_data.len() as f64 } else { f64::NAN };
    let w_data2: Vec<f64> = raw.iter().map(|r| r.w_1x2).filter(|x| x.is_finite()).collect();
    let w_1x2 = if w_data2.len()>0 { w_data2.iter().sum::<f64>()/w_data2.len() as f64 } else { f64::NAN };
    let w_data3: Vec<f64> = raw.iter().map(|r| r.w_2x2).filter(|x| x.is_finite()).collect();
    let w_2x2 = if w_data3.len()>0 { w_data3.iter().sum::<f64>()/w_data3.len() as f64 } else { f64::NAN };
    
    Meas { e:e_mean, e_err, v_abs:va, v_abs_err:va_err,
        v_stag:vs, v_stag_err:vs_err, binder:b, binder_err:b_err,
        tau_int_e:tau_e, gamma, l, beta:p.b,
        n_thermal, n_samples:cli.samples, n_spins:nspin(&p),
        wilson_1x1: w_1x1, wilson_1x2: w_1x2, wilson_2x2: w_2x2 }
}

// ============================================================
fn main() {
    let cli = Cli::parse();
    let gammas: Vec<f64> = if let Some(ref s) = cli.scan {
        s.split(',').map(|x| x.trim().parse().unwrap()).collect()
    } else { vec![cli.gamma] };
    let ls: Vec<usize> = if cli.fss { vec![4,6,8] } else { vec![cli.size] };
    
    if !cli.json {
        let header = format!("Ze QMC v1.2 | {} spins (i8) | PT={} | Wolff | τ_int{}",
            nspin(&Params{l:ls[0],lt:cli.lt,m:cli.trotter,jt:cli.jt,js:cli.js, jnnn:cli.jnnn,g:cli.gamma,h:cli.h,b:cli.beta}),
            cli.pt_replicas, if cli.auto_thermal {" | auto-thermal"} else {""});
        println!("{}", header);
        if cli.fss { println!("Binder crossing: L = {:?}", ls); }
        println!("{:>4} {:>6} {:>10} {:>10} {:>10} {:>10} {:>7} {:>5}",
                 "L","Γ","|v|","v_stag","E/N","Binder","τ_int","Фаза");
        println!("{}","─".repeat(72));
    }
    
    let mut all = vec![];
    for &l in &ls {
        for &g in &gammas {
            let m = run(&cli, g, l);
            let phase = if m.v_stag>0.3 {"АФМ"} else if m.v_abs<0.2 {"пара"} else {"крит"};
            // Wilson: area law (W_2x2≈W_1x1^4) → confinement; perimeter (W_2x2≈W_1x2^2) → deconfinement
            let w_diag = if m.wilson_1x1.is_finite() && m.wilson_2x2.is_finite() {
                let area_pred = m.wilson_1x1.powi(4);
                let perim_pred = m.wilson_1x2.powi(2);
                if (m.wilson_2x2-perim_pred).abs() < (m.wilson_2x2-area_pred).abs() {"deconf"} else {"conf"}
            } else {"?"};
            if cli.json { all.push(m.clone()); }
            else { println!("{:4} {:6.2} {:10.4} {:10.4} {:10.4} {:10.4} {:7.2} {:>5}  W:{:.3}/{:.3} {}",
                l,g,m.v_abs,m.v_stag,m.e,m.binder,m.tau_int_e,phase,m.wilson_1x1,m.wilson_2x2,w_diag); }
            all.push(m);
        }
    }
    
    // Binder crossing analysis
    if cli.fss && !cli.json && gammas.len()>=2 {
        println!("\n═══ Binder crossing analysis ═══");
        for &g in &gammas {
            let binders: Vec<(usize,f64)> = all.iter().filter(|m| (m.gamma-g).abs()<0.001)
                .map(|m| (m.l, m.binder)).collect();
            if binders.len()>=2 {
                print!("Γ={:.2}: ", g);
                for (l,b) in &binders { print!("L={}:{:.4} ", l, b); }
                let b_vals: Vec<f64> = binders.iter().map(|(_,b)| *b).collect();
                if b_vals.iter().all(|x| x.is_finite()) {
                    let b0 = b_vals[0]; let b_last = b_vals[b_vals.len()-1];
                    if (b0-0.66).abs()<0.1 && b_last<b0 { println!("→ crossing ✓ (ordered→disordered)"); }
                    else if b_last>b0 { println!("→ increasing (disordered→ordered)"); }
                    else { println!("→ flat"); }
                } else { println!(); }
            }
        }
    }
    
    if cli.json { println!("{}", serde_json::to_string_pretty(&all).unwrap()); }
}
