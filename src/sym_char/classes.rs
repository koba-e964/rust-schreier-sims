//! Conjugacy classes and related definitions.

use crate::perm::{all_permutations, get_cycle, Perm};
use num_bigint::BigInt;
use std::collections::HashSet;

pub trait ClassLike {
    fn to_repr(&self) -> Perm;
    fn get_class(p: &Perm) -> Self;
    /// Returns the size of the conjugacy class this represents.
    fn get_size(&self) -> BigInt;
}

impl ClassLike for Perm {
    fn to_repr(&self) -> Perm {
        self.clone()
    }
    fn get_class(p: &Perm) -> Self {
        p.clone()
    }
    fn get_size(&self) -> BigInt {
        1.into()
    }
}

/// A conjugacy class of some symmetric group S_n.
#[derive(Eq, PartialEq, Hash, Clone, Debug)]
pub struct ConjugacyClass {
    // part should be sorted in the descending order.
    part: Vec<usize>,
}

impl ClassLike for ConjugacyClass {
    fn to_repr(&self) -> Perm {
        let mut cycles = vec![];
        let mut cur = 0;
        for i in 0..self.part.len() {
            cycles.push((cur..cur + self.part[i]).collect());
            cur += self.part[i];
        }
        get_cycle(cur, &cycles)
    }
    fn get_class(p: &Perm) -> Self {
        let n = p.size();
        let mut visited = vec![false; n];
        let mut part = vec![];
        for i in 0..n {
            if visited[i] {
                continue;
            }
            let mut cur = p[i];
            visited[i] = true;
            let mut length = 1;
            while cur != i {
                visited[cur] = true;
                cur = p[cur];
                length += 1;
            }
            part.push(length);
        }
        part.sort();
        part.reverse();
        Self { part }
    }
    fn get_size(&self) -> BigInt {
        let mut prod = BigInt::from(1);
        let n: usize = self.part.iter().sum();
        for i in 1..n + 1 {
            prod *= i;
        }
        let mut pre = 0;
        let mut count = 0;
        for &v in &self.part {
            prod /= v;
            if pre != v {
                count = 0;
            }
            pre = v;
            count += 1;
            prod /= count;
        }
        prod
    }
}

/// Returns all conjugacy classes of S_n.
/// It is guaranteed that returned[0] contains e.
pub fn get_all_concjugacy_classes(n: usize) -> Vec<ConjugacyClass> {
    // TODO: just enumerate all partitions of n, instead of enumerating all n! permutations
    let perms = all_permutations(n);
    let mut conjs = vec![];
    let mut seen = HashSet::new();
    for p in perms {
        let conj = ConjugacyClass::get_class(&p);
        if seen.contains(&conj) {
            continue;
        }
        seen.insert(conj.clone());
        conjs.push(conj);
    }
    conjs
}
