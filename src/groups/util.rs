use crate::perm::Perm;

/// Given a list of cycles, this function returns its composition.
pub fn get_cycle(n: usize, a: &[Vec<usize>]) -> Perm {
    let mut e = Perm::e(n);
    for a in a {
        let mut t: Vec<_> = (0..n).collect();
        for i in 0..a.len() {
            t[a[i]] = a[(i + 1) % a.len()];
        }
        e = e.compose(&Perm::new(t));
    }
    e
}
