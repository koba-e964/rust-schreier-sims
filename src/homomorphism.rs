use crate::perm::Perm;
use crate::schreier::incrementally_build_bsgs;

/// Checks if a given mapping (x[i] |-> y[i]) extends to a group homomorphism
/// <x> -> <y>.
pub fn is_homomorphism(n: usize, m: usize, x: &[Perm], y: &[Perm]) -> bool {
    assert_eq!(x.len(), y.len());
    let mut rnd = rand::thread_rng();

    let mut xy = vec![Perm::e(1); x.len()];
    // Make a concatenation of x[i] and y[i] for each i.
    for i in 0..x.len() {
        xy[i] = x[i].concat(&y[i]);
    }
    let (beta_transversals, _) = incrementally_build_bsgs(n, &[], x, &mut rnd);
    let beta: Vec<_> = beta_transversals.iter().map(|&(b, _)| b).collect();
    eprintln!("beta = {:?}", beta);
    // Is beta also a BSGS of <xy>?
    let (beta_transversals, _) = incrementally_build_bsgs(n + m, &beta, &xy, &mut rnd);
    beta_transversals.len() == beta.len()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::groups::rubik;

    #[test]
    fn is_homomorphism_test_a2_a3() {
        let x = vec![Perm::new(vec![1, 0])];
        let y = vec![Perm::new(vec![1, 2, 0])];
        assert!(!is_homomorphism(2, 3, &x, &y));
    }

    #[test]
    fn is_homomorphism_test_d8_s2() {
        let x = vec![Perm::new(vec![1, 2, 3, 0]), Perm::new(vec![2, 1, 0, 3])];
        let y = vec![Perm::new(vec![1, 0]), Perm::e(2)];
        assert!(is_homomorphism(4, 2, &x, &y));
    }
    #[test]
    fn is_homomorphism_test_s6_outer_automorphism() {
        // _The_ outer automorphism S6 -> S6
        // described in http://www2u.biglobe.ne.jp/~nuida/m/doc/OutS6.pdf.
        // This test verifies that the map given in the article is indeed
        // a homomorphism.
        let transposition = |x, y| {
            let mut v: Vec<usize> = (0..6).collect();
            v.swap(x, y);
            Perm::new(v)
        };
        let x = vec![
            transposition(0, 1),
            transposition(1, 2),
            transposition(2, 3),
            transposition(3, 4),
            transposition(4, 5),
        ];
        let y = vec![
            Perm::new(vec![1, 0, 3, 2, 5, 4]),
            Perm::new(vec![5, 3, 4, 1, 2, 0]),
            Perm::new(vec![3, 2, 1, 0, 5, 4]),
            Perm::new(vec![5, 4, 3, 2, 1, 0]),
            Perm::new(vec![2, 3, 0, 1, 5, 4]),
        ];
        assert!(is_homomorphism(6, 6, &x, &y));
    }
    #[test]
    fn is_homomorphism_test_rubik_corner_cubes() {
        // Rubik's Cube group
        let (n, gen) = rubik::generators();
        assert_eq!(n, 48);
        // Corner cubes are labeled with even numbers
        let evengen: Vec<_> = gen
            .iter()
            .map(|g| {
                let m = 24;
                let mut v = vec![0; m];
                for i in 0..m {
                    v[i] = g[2 * i] / 2;
                }
                Perm::new(v)
            })
            .collect();
        is_homomorphism(n, n / 2, &gen, &evengen);
    }
    #[test]
    fn is_homomorphism_test_rubik_edge_cubes() {
        // Rubik's Cube group
        let (n, gen) = rubik::generators();
        assert_eq!(n, 48);
        // Edge cubes are labeled with odd numbers
        let oddgen: Vec<_> = gen
            .iter()
            .map(|g| {
                let m = 24;
                let mut v = vec![0; m];
                for i in 0..m {
                    v[i] = g[2 * i + 1] / 2;
                }
                Perm::new(v)
            })
            .collect();
        is_homomorphism(n, n / 2, &gen, &oddgen);
    }
}
