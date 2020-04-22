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

    pub fn sgn(&self) -> i8 {
        // inversion count, O(n^2)
        let n = self.size();
        let mut ans = 1;
        for i in 0..n {
            for j in 0..i {
                if self[j] > self[i] {
                    ans = -ans;
                }
            }
        }
        ans
    }
}

impl std::ops::Index<usize> for Perm {
    type Output = usize;

    fn index(&self, index: usize) -> &usize {
        let Perm(x) = &self;
        &x[index]
    }
}

pub fn all_permutations(n: usize) -> Vec<Perm> {
    let mut ans = vec![];

    let mut v: Vec<usize> = (0..n).collect();

    loop {
        ans.push(Perm::new(v.clone()));
        let mut tail_dec: usize = 1;
        while tail_dec < n {
            if v[n - tail_dec - 1] > v[n - tail_dec] {
                tail_dec += 1;
            } else {
                break;
            }
        }
        // v[n - tail_dec .. n] is strictly decreasing
        if tail_dec < n {
            let x = n - tail_dec - 1;
            let mut y = n;
            {
                let pivot = &v[x];
                for i in (n - tail_dec..n).rev() {
                    if v[i] > *pivot {
                        y = i;
                        break;
                    }
                }
                assert!(y < n);
            }
            v.swap(x, y);
        }
        v[n - tail_dec..n].reverse();
        if tail_dec >= n {
            break;
        }
    }
    ans
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
    fn all_permutations_test() {
        let all = all_permutations(3);
        assert_eq!(all.len(), 6);
    }
}
