use crate::perm::Perm;
use crate::transversal::{get_transversal, orbit_transversal_stabilizer, Transversal};
use rand::Rng;

// Reference: https://blogs.cs.st-andrews.ac.uk/codima/files/2015/11/CoDiMa2015_Holt.pdf

/// The returned vaule ([u0, ...], h) must satisfy g = u0 * ... * h.
pub fn strip(g: &Perm, beta_transversals: &[(usize, Transversal)]) -> (Vec<Perm>, Perm) {
    let mut h = g.clone();
    let mut us = vec![];
    for &(beta, ref transversal) in beta_transversals {
        let moved_to = h[beta];
        match transversal[moved_to] {
            // If repr is dummy, that is, moved_to is not in the orbit beta^H
            None => break,
            Some(ref repr) => {
                h = h.compose(&repr.inv());
                us.push(repr.clone());
            }
        }
    }
    (us, h)
}

pub fn schreier_sims(
    n: usize,
    beta_transversals: &[(usize, Transversal)],
    s: &[Perm],
) -> Result<(), (Vec<Perm>, Perm)> {
    if beta_transversals.is_empty() {
        if s.is_empty() {
            return Ok(());
        }
        return Err((vec![], s[0].clone()));
    }
    let mut intersection = Vec::new();
    // The first fixed point
    let beta0 = beta_transversals[0].0;
    for s in s {
        if s[beta0] == beta0 {
            intersection.push(s.clone());
        }
    }
    schreier_sims(n, &beta_transversals[1..], &intersection)?;
    let (_, y) = orbit_transversal_stabilizer(n, s, beta0);
    for y in y {
        let (us, rest) = strip(&y, &beta_transversals[1..]);
        if rest != Perm::e(n) {
            return Err((us, rest));
        }
    }
    Ok(())
}

/// Returns B and S built.
pub fn incrementally_build_bsgs(
    n: usize,
    initial_beta: &[usize],
    initial_s: &[Perm],
    mut rnd: impl Rng,
) -> (Vec<(usize, Transversal)>, Vec<Perm>) {
    let mut beta_transversals = vec![];
    let mut s = initial_s.to_vec();
    let mut used = vec![false; n];
    let dummy_transversal = vec![]; // dummy transversals
    for &beta in initial_beta {
        beta_transversals.push((beta, dummy_transversal.clone()));
    }
    loop {
        // preliminary result
        // TODO reuse transversals
        {
            let mut cur_s = s.clone();
            for &mut (beta, ref mut transversal_ref) in &mut beta_transversals {
                let (orbit_transversal, _) = orbit_transversal_stabilizer(n, &cur_s, beta);
                let transversal = get_transversal(n, orbit_transversal);
                *transversal_ref = transversal;
                cur_s.retain(|perm| perm[beta] == beta);
                used[beta] = true;
            }
        }

        // Incrementally computes Y and check if it's okay.
        match schreier_sims(n, &beta_transversals, &s) {
            Ok(()) => break,
            Err((_, h)) => {
                s.push(h.clone());
                // Are there any points that are not stabilized by h
                // and in beta?
                let mut moved = vec![];
                for i in 0..n {
                    if h[i] != i {
                        if used[i] {
                            moved.clear();
                            break;
                        }
                        moved.push(i);
                    }
                }
                if !moved.is_empty() {
                    // All points that are not stabilized by h are not in beta.
                    // randomly pick one of them
                    let point = moved[rnd.gen_range(0..moved.len())];
                    beta_transversals.push((point, dummy_transversal.clone()));
                }
            }
        }
    }
    (beta_transversals, s)
}

pub fn order(n: usize, gen: &[Perm]) -> num_bigint::BigInt {
    let mut order = 1.into();
    let mut rnd = rand::thread_rng();
    let (beta_transversals, _) = incrementally_build_bsgs(n, &[], gen, &mut rnd);
    for (_, transversal) in beta_transversals {
        let mut u = 0;
        for transversal in transversal {
            if transversal.is_some() {
                u += 1;
            }
        }
        order *= u;
    }
    order
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::groups::{mathieu12, rubik};

    #[test]
    fn schreier_sims_test() {
        // G= <(0 1 2), (2 3 4)>
        let n = 5;
        let gen = vec![
            Perm::new(vec![1, 2, 0, 3, 4]),
            Perm::new(vec![0, 1, 3, 4, 2]),
        ];
        let beta = vec![0, 2];
        let (orbit_transversal0, subgen) = orbit_transversal_stabilizer(n, &gen, beta[0]);
        assert_eq!(orbit_transversal0.len(), 5);
        let transversal0 = get_transversal(n, orbit_transversal0);
        let (orbit_transversal1, _) = orbit_transversal_stabilizer(n, &subgen, beta[1]);
        assert_eq!(orbit_transversal1.len(), 4);
        let transversal1 = get_transversal(n, orbit_transversal1);
        // orbit 1^{G^{(1)}} is {1, 2, 3, 4}.
        assert_eq!(transversal1[0], None);

        let beta_transversals = vec![(beta[0], transversal0), (beta[1], transversal1)];
        let ans = schreier_sims(5, &beta_transversals, &gen);
        // (beta, gen) is not a BSGS:
        // there is an element of G that needs appending to gen.
        assert!(ans.is_err());
        let (_, h) = ans.unwrap_err();
        // Asserts h stabilizes 0 and 2.
        for &point in &beta {
            assert_eq!(h[point], point);
        }
    }
    #[test]
    fn incrementally_build_bsgs_test() {
        // G = <(0 1 2), (2 3 4)>
        let n = 5;
        let gen = vec![
            Perm::new(vec![1, 2, 0, 3, 4]),
            Perm::new(vec![0, 1, 3, 4, 2]),
        ];
        let beta = vec![0, 2];
        let mut rnd = rand::thread_rng();
        let (beta_transversals, _) = incrementally_build_bsgs(n, &beta, &gen, &mut rnd);
        // This result is hard to test,
        // but it can be tested by computing |G| = \Prod |U_i|.
        // |G| should be 60. In fact, G = A_5.
        let mut order = 1;
        for (_, transversal) in beta_transversals {
            let mut u = 0;
            for i in 0..n {
                if let Some(_) = transversal[i] {
                    u += 1;
                }
            }
            order *= u;
        }
        assert_eq!(order, 60);
    }
    #[test]
    fn order_test_0() {
        // G = <(0 1 2), (2 3 4)> = A_5, |G| = 60
        let n = 5;
        let gen = vec![
            Perm::new(vec![1, 2, 0, 3, 4]),
            Perm::new(vec![0, 1, 3, 4, 2]),
        ];
        assert_eq!(order(n, &gen), 60.into());
    }
    #[test]
    fn order_test_1() {
        // G = <(0 1 2 3), (0 2)> ~= D_8, |G| = 8
        let n = 4;
        let gen = vec![Perm::new(vec![1, 2, 3, 0]), Perm::new(vec![2, 1, 0, 3])];
        assert_eq!(order(n, &gen), 8.into());
    }
    #[test]
    fn order_test_2() {
        // Star, G = <(0 n-1), (1 n-1), ...> = S_n, |G| = n!
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
        let mut factorial: num_bigint::BigInt = 1.into();
        for i in 1..n + 1 {
            factorial *= i;
        }
        assert_eq!(order(n, &gen), factorial);
    }
    #[test]
    fn order_test_3() {
        use std::str::FromStr;
        // Rubik's Cube group
        let (n, gen) = rubik::generators();
        assert_eq!(
            order(n, &gen),
            num_bigint::BigInt::from_str("43252003274489856000").unwrap()
        );
    }
    #[test]
    fn order_test_4() {
        // The Mathieu group M12
        let (n, gen) = mathieu12::generators();
        assert_eq!(order(n, &gen), 95040.into());
    }
}
