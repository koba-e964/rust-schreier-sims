use rand::Rng;
use crate::perm::Perm;
use crate::transversal::{Transversal, orbit_transversal_stabilizer};

// Reference: https://blogs.cs.st-andrews.ac.uk/codima/files/2015/11/CoDiMa2015_Holt.pdf

pub fn strip(
    g: &Perm,
    beta_transversals: &[(usize, Transversal)]
) -> (Vec<Perm>, Perm) {
    let mut h = g.clone();
    let mut us = vec![];
    for i in 0..beta_transversals.len() {
        let (beta, ref transversal) = beta_transversals[i];
        let moved_to = h[beta];
        let repr = &transversal[moved_to];
        // If repr is dummy, that is, moved_to is not in beta^H
        // TODO refactor
        if repr.size() == 0 {
            return (us, h);
        }
        h = h.compose(&repr.inv());
        us.push(repr.clone());
    }
    (us, h)
}

pub fn schreier_sims(
    n: usize,
    beta_transversals: &[(usize, Transversal)],
    s: &[Perm],
) -> Result<(), (Vec<Perm>, Perm)> {
    if beta_transversals.len() == 0 {
        if s.len() == 0 {
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
    let (_orbit, y) = orbit_transversal_stabilizer(n, s, beta0);
    for y in y {
        let (u, h) = strip(&y, &beta_transversals[1..]);
        if h != Perm::e(n) {
            return Err((u, h));
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
    for i in 0..initial_beta.len() {
        beta_transversals.push((initial_beta[i], dummy_transversal.clone()));
    }
    loop {
        // preliminary result
        // TODO reuse transversals
        {
            let mut cur_s = s.clone();
            for i in 0..beta_transversals.len() {
                let beta = beta_transversals[i].0;
                let (orbit_transversal, _) = orbit_transversal_stabilizer(n, &cur_s, beta);
                let mut transversal = vec![Perm::e(0); n];
                for (point, trans) in orbit_transversal {
                    transversal[point] = trans;
                }
                beta_transversals[i].1 = transversal;
                cur_s = cur_s
                    .into_iter()
                    .filter(|perm| perm[beta] == beta)
                    .collect();
                used[beta] = true;
            }
        }

        // Incrementally computes Y and check if it's okay.
        match schreier_sims(n, &beta_transversals, &s) {
            Ok(()) =>
                break,
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
                if moved.len() > 0 {
                    // All points that are not stabilized by h are not in beta.
                    // randomly pick one of them
                    let point = moved[rnd.gen_range(0, moved.len())];
                    beta_transversals.push((point, dummy_transversal.clone()));
                }
            }
        }
    }
    (beta_transversals, s)
}

#[cfg(test)]
mod tests {
    use super::*;

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
        let mut transversal0 = vec![Perm::e(0); n];
        for (point, trans) in orbit_transversal0 {
            transversal0[point] = trans;
        }
        let (orbit_transversal1, _) = orbit_transversal_stabilizer(n, &subgen, beta[1]);
        assert_eq!(orbit_transversal1.len(), 4);
        let mut transversal1 = vec![Perm::e(0); n];
        for (point, trans) in orbit_transversal1 {
            transversal1[point] = trans;
        }
        // orbit 1^{G^{(1)}} is {1, 2, 3, 4}.
        assert_eq!(transversal1[0], Perm::e(0));

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
                if transversal[i].size() != 0 {
                    u += 1;
                }
            }
            order *= u;
        }
        assert_eq!(order, 60);
    }
}
