#[macro_use]
extern crate criterion;

use criterion::Criterion;
use rust_schreier::perm::Perm;
use rust_schreier::schreier::order;

fn star_benchmark_10(c: &mut Criterion) {
    // Star, G = <(0 n-1), (1 n-1), ...> = S_n, |G| = n!
    c.bench_function("star 10", |b| {
        b.iter(|| {
            let n = 10;
            let mut gen = vec![];
            for i in 0..n - 1 {
                let mut p = vec![0; n];
                for j in 0..n {
                    p[j] = j;
                }
                p.swap(i, n - 1);
                gen.push(Perm::new(p));
            }
            order(n, &gen)
        })
    });
}
fn star_benchmark_20(c: &mut Criterion) {
    // Star, G = <(0 n-1), (1 n-1), ...> = S_n, |G| = n!
    c.bench_function("star 20", |b| {
        b.iter(|| {
            let n = 20;
            let mut gen = vec![];
            for i in 0..n - 1 {
                let mut p = vec![0; n];
                for j in 0..n {
                    p[j] = j;
                }
                p.swap(i, n - 1);
                gen.push(Perm::new(p));
            }
            order(n, &gen)
        })
    });
}
fn star_benchmark_30(c: &mut Criterion) {
    // Star, G = <(0 n-1), (1 n-1), ...> = S_n, |G| = n!
    c.bench_function("star 30", |b| {
        b.iter(|| {
            let n = 30;
            let mut gen = vec![];
            for i in 0..n - 1 {
                let mut p = vec![0; n];
                for j in 0..n {
                    p[j] = j;
                }
                p.swap(i, n - 1);
                gen.push(Perm::new(p));
            }
            order(n, &gen)
        })
    });
}
criterion_group!(
    benches,
    star_benchmark_10,
    star_benchmark_20,
    star_benchmark_30
);
criterion_main!(benches);
