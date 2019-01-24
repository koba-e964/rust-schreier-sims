use std::collections::VecDeque;

// Reference: https://blogs.cs.st-andrews.ac.uk/codima/files/2015/11/CoDiMa2015_Holt.pdf
#[derive(Clone, Eq, PartialEq, Debug, PartialOrd, Ord)]
struct Perm(Vec<usize>);
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
    fn pow(&self, k: i64) -> Self {
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

/// gen: generators, v: stabilized point
///
/// This function returns a generator set of the stabilizer group G_v = Stab_G(v).
fn schreier_theorem(n: usize, gen: &[Perm], v: usize) -> (usize, Vec<Perm>) {
    let mut stabilizer_gen = Vec::new();
    // Calculates a variant of the Schreier vector
    // performing breadth-first search.
    let mut que = VecDeque::new();
    que.push_back((v, Perm::e(n)));
    let mut table: Vec<Option<Perm>> = vec![None; n];
    let mut orbit = 0;
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
            table[w] = Some(p);
            orbit += 1;
        }
    }
    // Remove duplicate generators
    stabilizer_gen.sort();
    stabilizer_gen.dedup();

    // Returns ((G : G_v), a generator set of G_v)
    (orbit, stabilizer_gen)
}

fn main() {
    let gen = vec![Perm::new(vec![1, 2, 3, 0]), Perm::new(vec![1, 0, 2, 3])];
    let (h1, stab1) = schreier_theorem(4, &gen, 0);
    println!("{:?}", (h1, &stab1));
    let (h2, stab2) = schreier_theorem(4, &stab1, 1);
    println!("{:?}", (h2, &stab2));
    let (h3, stab3) = schreier_theorem(4, &stab2, 2);
    println!("{:?}", (h3, &stab3));
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
    fn schreier_theorem_test() {
        let gen = vec![Perm::new(vec![1, 2, 3, 0]), Perm::new(vec![1, 0, 2, 3])];
        let (h1, stab1) = schreier_theorem(4, &gen, 0);
        // |0^G| = 4, H^{(1)} = G_0 = <(2 3), (1 2)>
        assert_eq!(h1, 4);
        let (h2, stab2) = schreier_theorem(4, &stab1, 1);
        // |1^{H^{(1)}}| = 3, H^{(2)} = H^{(1)}_1 = <(2 3)>
        assert_eq!(h2, 3);
        let (h3, stab3) = schreier_theorem(4, &stab2, 2);
        // |2^{H^{(2)}}| = 2, H^{(3)} = H^{(2)}_2 = {e}
        assert_eq!(h3, 2);
        assert_eq!(stab3, Vec::new());
    }
}
