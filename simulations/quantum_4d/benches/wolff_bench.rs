use criterion::{black_box, Criterion, criterion_group, criterion_main};
use ze_qmc_4d::*;
use rand_xoshiro::Xoshiro256PlusPlus;
use rand::SeedableRng;

fn bench_wolff(c: &mut Criterion) {
    let p = Params { l:4, lt:6, m:16, jt:1.0, js:0.0, jnnn:0.0, g:1.0, h:0.0, b:10.0 };
    let cpl = TC::new(&p);
    let mut z = init_staggered(&p);
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(42);
    
    c.bench_function("wolff L=4", |b| {
        b.iter(|| { wolff(black_box(&mut z), black_box(&p), black_box(&cpl), black_box(&mut rng)); })
    });
}

fn bench_energy(c: &mut Criterion) {
    let p = Params { l:4, lt:6, m:16, jt:1.0, js:0.0, jnnn:0.0, g:1.0, h:0.0, b:10.0 };
    let cpl = TC::new(&p);
    let z = init_staggered(&p);
    
    c.bench_function("energy L=4", |b| {
        b.iter(|| { energy_config(black_box(&z), black_box(&p), black_box(&cpl)); })
    });
}

fn bench_measure(c: &mut Criterion) {
    let p = Params { l:4, lt:6, m:16, jt:1.0, js:0.0, jnnn:0.0, g:1.0, h:0.0, b:10.0 };
    let cpl = TC::new(&p);
    let z = init_staggered(&p);
    
    c.bench_function("measure L=4", |b| {
        b.iter(|| { measure_one(black_box(&z), black_box(&p), black_box(&cpl)); })
    });
}

criterion_group!(benches, bench_wolff, bench_energy, bench_measure);
criterion_main!(benches);
