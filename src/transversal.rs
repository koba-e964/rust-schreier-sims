use std::collections::VecDeque;
use crate::perm::Perm;

// Reference: https://blogs.cs.st-andrews.ac.uk/codima/files/2015/11/CoDiMa2015_Holt.pdf

// A collection of pairs (x, alpha) s.t. v^alpha = x.
pub type OrbitTransversal = Vec<(usize, Perm)>;

// An oracle of type usize -> Perm.
pub type Transversal = Vec<Perm>;

/// gen: generators, v: stabilized point
///
/// This function returns a generator set of the stabilizer group G_v = Stab_G(v).
pub fn orbit_transversal_stabilizer(n: usize, gen: &[Perm], v: usize) -> (OrbitTransversal, Vec<Perm>) {
    let mut stabilizer_gen = Vec::new();
    // Calculates a variant of the Schreier vector
    // performing breadth-first search.
    let mut que = VecDeque::new();
    que.push_back((v, Perm::e(n)));
    let mut table: Vec<Option<Perm>> = vec![None; n];
    let mut orbit_transversal = vec![];
    while let Some((w, p)) = que.pop_front() {
        if let Some(ref q) = table[w] {
            // r = qp^{-1} stabilizes v.
            let r = q.compose(&p.inv());
            if r != Perm::e(n) {
                stabilizer_gen.push(r);
            }
        } else {
            for x in gen {
                // move to w^x (x applied to w).
                que.push_back((x[w], p.compose(x)));
            }
            orbit_transversal.push((w, p.clone()));
            table[w] = Some(p);
        }
    }
    // Remove duplicate generators
    stabilizer_gen.sort();
    stabilizer_gen.dedup();

    // Returns ((G : G_v), a generator set of G_v)
    (orbit_transversal, stabilizer_gen)
}

pub fn strip(g: &Perm, beta_transversals: &[(usize, Transversal)]) -> (Vec<Perm>, Perm) {
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

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn orbit_transversal_stabilizer_test() {
        let gen = vec![Perm::new(vec![1, 2, 3, 0]), Perm::new(vec![1, 0, 2, 3])];
        let (h1, stab1) = orbit_transversal_stabilizer(4, &gen, 0);
        // |0^G| = 4, H^{(1)} = G_0 = <(2 3), (1 2)>
        assert_eq!(h1.len(), 4);
        let (h2, stab2) = orbit_transversal_stabilizer(4, &stab1, 1);
        // |1^{H^{(1)}}| = 3, H^{(2)} = H^{(1)}_1 = <(2 3)>
        assert_eq!(h2.len(), 3);
        let (h3, stab3) = orbit_transversal_stabilizer(4, &stab2, 2);
        // |2^{H^{(2)}}| = 2, H^{(3)} = H^{(2)}_2 = {e}
        assert_eq!(h3.len(), 2);
        assert_eq!(stab3, Vec::new());
    }
}
