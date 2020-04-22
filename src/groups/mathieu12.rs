// Utility functions for the Mathieu group M12.
use crate::groups::util::get_cycle;
use crate::perm::Perm;

/// Reference: http://brauer.maths.qmul.ac.uk/Atlas/v3/permrep/M12G1-p12aB0
pub fn generators() -> (usize, Vec<Perm>) {
    // In GAP, the following two generators are defined by:
    // b11 := (1,4)(3,10)(5,11)(6,12);
    // b21 := (1,8,9)(2,3,4)(5,12,11)(6,10,7);
    let b11 = get_cycle(12, &[vec![0, 3], vec![2, 9], vec![4, 10], vec![5, 11]]);
    let b21 = get_cycle(
        12,
        &[vec![0, 7, 8], vec![1, 2, 3], vec![4, 11, 10], vec![5, 9, 6]],
    );
    (12, vec![b11, b21])
}
