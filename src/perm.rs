#[derive(Clone, Eq, PartialEq, Debug, PartialOrd, Ord)]
pub struct Perm(Vec<usize>);
impl Perm {
    pub fn new(perm: Vec<usize>) -> Perm {
        let n = perm.len();
        let mut appear = vec![false; n];
        for &p in &perm {
            assert!(!appear[p]);
            appear[p] = true;
        }
        Perm(perm)
    }
    pub fn size(&self) -> usize {
        let Perm(x) = self;
        x.len()
    }
    pub fn compose(&self, Perm(other): &Self) -> Self {
        let Perm(me) = self;
        let n = self.size();
        assert_eq!(n, other.len());
        let mut ans = vec![0; n];
        for i in 0..n {
            ans[i] = other[me[i]];
        }
        Self::new(ans)
    }
    pub fn e(n: usize) -> Self {
        Self::new((0..n).collect())
    }
    pub fn inv(&self) -> Self {
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
    /// Returns a new permutation of self.size() + a.size().
    pub fn concat(&self, a: &Perm) -> Perm {
        let n = self.size();
        let m = a.size();
        let mut v = vec![0; n + m];

        for i in 0..n {
            v[i] = self[i];
        }
        for i in 0..m {
            v[n + i] = n + a[i];
        }
        Perm::new(v)
    }
}

impl std::ops::Index<usize> for Perm {
    type Output = usize;

    fn index(&self, index: usize) -> &usize {
        let Perm(x) = &self;
        &x[index]
    }
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
}
