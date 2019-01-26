use rust_schreier::perm::Perm;
use rust_schreier::schreier::order;

fn main() {
    // Star, G = <(0 n-1), (1 n-1), ...> = S_n, |G| = n!
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
    eprintln!("order = {}", order(n, &gen));
}
