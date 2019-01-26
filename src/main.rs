use std::collections::VecDeque;

// Reference: https://blogs.cs.st-andrews.ac.uk/codima/files/2015/11/CoDiMa2015_Holt.pdf
#[derive(Clone, Eq, PartialEq, Debug, PartialOrd, Ord)]
pub struct Perm(Vec<usize>);
impl Perm {
    fn new(perm: Vec<usize>) -> Perm {
        let n = perm.len();
        let mut appear = vec![false; n];
        for &p in &perm {
            assert!(!appear[p]);
            appear[p] = true;
        }
        Perm(perm)
    }
    fn size(&self) -> usize {
        let Perm(x) = self;
        x.len()
    }
    fn compose(&self, Perm(other): &Self) -> Self {
        let Perm(me) = self;
        let n = self.size();
        assert_eq!(n, other.len());
        let mut ans = vec![0; n];
        for i in 0..n {
            ans[i] = other[me[i]];
        }
        Self::new(ans)
    }
    fn e(n: usize) -> Self {
        Self::new((0..n).collect())
    }
    fn inv(&self) -> Self {
        let Perm(me) = self;
        let n = self.size();
        let mut ans = vec![0; n];
        for i in 0..n {
            ans[me[i]] = i;
        }
        Self::new(ans)
    }
    pub fn pow(&self, k: i64) -> Self {
        let n = self.size();
        let mut sum = Self::e(n);
        if k == 0 {
            return sum;
        }
        let mut cur = if k < 0 { self.inv() } else { self.clone() };
        let mut k = k.abs();
        while k > 0 {
            if (k & 1) == 1 {
                sum = sum.compose(&cur);
            }
            cur = cur.compose(&cur);
            k >>= 1;
        }
        sum
    }
}

impl std::ops::Index<usize> for Perm {
    type Output = usize;

    fn index(&self, index: usize) -> &usize {
        let Perm(x) = &self;
        &x[index]
    }
}

// A collection of pairs (x, alpha) s.t. v^alpha = x.
type OrbitTransversal = Vec<(usize, Perm)>;

// An oracle of type usize -> Perm.
type Transversal = Vec<Perm>;

/// gen: generators, v: stabilized point
///
/// This function returns a generator set of the stabilizer group G_v = Stab_G(v).
fn orbit_transversal_stabilizer(n: usize, gen: &[Perm], v: usize)
                                -> (OrbitTransversal, Vec<Perm>) {
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

fn strip(g: &Perm, beta_transversals: &[(usize, Transversal)]) -> (Vec<Perm>, Perm) {
    let mut h = g.clone();
    let mut us = vec![];
    for i in 0 .. beta_transversals.len() {
        let (beta, ref transversal) = beta_transversals[i];
        let moved_to = h[beta];
        eprintln!("trans = {:?}, {} -> {}", transversal, beta, moved_to);
        let repr = &transversal[moved_to];
        h = h.compose(&repr.inv());
        us.push(repr.clone());
    }
    (us, h)
}

fn schreier_sims(n: usize, beta_transversals: &[(usize, Transversal)],
                 s: &[Perm]) -> Result<(), (Vec<Perm>, Perm)> {
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

fn main() {
    // G= <(0 1 2), (2 3 4)>
    let n = 5;
    let gen = vec![Perm::new(vec![1, 2, 0, 3, 4]),
                   Perm::new(vec![0, 1, 3, 4, 2])];
    let beta = vec![0, 2];
    let (orbit_transversal0, subgen) =
        orbit_transversal_stabilizer(n, &gen, beta[0]);
    let mut transversal0 = vec![Perm::e(0); n];
    for (point, trans) in orbit_transversal0 {
        transversal0[point] = trans;
    }
    eprintln!("transversal0 = {:?}", transversal0);
    let (orbit_transversal1, _) =
        orbit_transversal_stabilizer(n, &subgen, beta[1]);
    let mut transversal1 = vec![Perm::e(0); n];
    for (point, trans) in orbit_transversal1 {
        transversal1[point] = trans;
    }
    eprintln!("transversal1 = {:?}", transversal1);

    let beta_transversals = vec![(beta[0], transversal0), (beta[1], transversal1)];
    let ans = schreier_sims(5, &beta_transversals, &gen);
    eprintln!("ans = {:?}", ans);
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn perm_pow_test() {
        let p = Perm::new(vec![1, 2, 3, 4, 0]);
        // p^8 = [3, 4, 0, 1, 2]
        assert_eq!(p.pow(8), Perm::new(vec![3, 4, 0, 1, 2]));
        // p^5 = e
        assert_eq!(p.pow(5), Perm::e(5));
        // p^{-5} = e
        assert_eq!(p.pow(-5), Perm::e(5));
        // p^{-3} = p^2 = pp
        assert_eq!(p.pow(-3), p.compose(&p));
    }
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
    #[test]
    fn schreier_sims_test() {
        // G= <(0 1 2), (2 3 4)>
        let n = 5;
        let gen = vec![Perm::new(vec![1, 2, 0, 3, 4]),
                       Perm::new(vec![0, 1, 3, 4, 2])];
        let beta = vec![0, 2];
        let (orbit_transversal0, subgen) =
            orbit_transversal_stabilizer(n, &gen, beta[0]);
        assert_eq!(orbit_transversal0.len(), 5);
        let mut transversal0 = vec![Perm::e(0); n];
        for (point, trans) in orbit_transversal0 {
            transversal0[point] = trans;
        }
        let (orbit_transversal1, _) =
            orbit_transversal_stabilizer(n, &subgen, beta[1]);
        assert_eq!(orbit_transversal1.len(), 4);
        let mut transversal1 = vec![Perm::e(0); n];
        for (point, trans) in orbit_transversal1 {
            transversal1[point] = trans;
        }
        // orbit 1^{G^{(1)}} is {1, 2, 3, 4}.
        assert_eq!(transversal1[0], Perm::e(0));
        
        let beta_transversals =
            vec![(beta[0], transversal0), (beta[1], transversal1)];
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
}
