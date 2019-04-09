use crate::perm::Perm;
use crate::schreier::order;

/// Checks if a given mapping (x[i] |-> y[i]) extends to a group homomorphism
/// <x> -> <y>.
pub fn is_homomorphism(n: usize, m: usize, x: &[Perm], y: &[Perm]) -> bool {
    assert_eq!(x.len(), y.len());
    let g = order(n, x);
    eprintln!("g = {}", g);
    let mut xy = vec![Perm::e(1); x.len()];
    for i in 0..x.len() {
        // Make a concatenation of x[i] and y[i].
        let mut p = vec![0; n + m];
        for j in 0..n {
            p[j] = x[i][j];
        }
        for j in 0..m {
            p[n + j] = y[i][j] + n;
        }
        xy[i] = Perm::new(p);
    }
    let gh = order(n + m, &xy);
    eprintln!("gh = {}", gh);
    g == gh
}

#[cfg(test)]
mod tests {
    use super::*;

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
}
