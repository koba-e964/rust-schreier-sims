use rust_schreier::perm::Perm;
use rust_schreier::schreier::incrementally_build_bsgs;


fn main() {
    // G = <(0 1 2), (2 3 4)>
    let n = 5;
    let gen = vec![
        Perm::new(vec![1, 2, 0, 3, 4]),
        Perm::new(vec![0, 1, 3, 4, 2]),
    ];
    let beta = vec![0, 2];
    let mut rnd = rand::thread_rng();
    let (beta_transversals, s) = incrementally_build_bsgs(n, &beta, &gen, &mut rnd);
    for (beta, transversal) in beta_transversals {
        eprintln!("beta: {}, transversal = {:?}", beta, transversal);
    }
    eprintln!("s = {:?}", s);
}
