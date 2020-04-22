use char_entry::{CharEntry, ClassLike};
use num_bigint::BigInt;

mod char_entry;

pub struct CharTableSummary {
    pub table: Vec<CharEntry>,
    pub remaining: BigInt,
}

fn find_new_entry<C: ClassLike>(g: &[C], ans: &[CharEntry]) -> Option<CharEntry> {
    let n = ans.len();

    for i in 0..n {
        for j in i..n {
            let mut tensor = ans[i].tensor_prod(&ans[j]);
            for j in 0..n {
                let inp = tensor.inner_prod(&ans[j], g);
                tensor -= &ans[j].power(inp);
            }
            let norm = tensor.inner_prod(&tensor, g);
            eprintln!(
                "found: ({}) tensor ({}) gives a rep of dim {}, norm {}",
                i, j, tensor.table[0], norm
            );
            if norm == 1.into() {
                // A new irrep was found.
                return Some(tensor);
            }
        }
    }
    None
}

pub fn get_char_table<C: ClassLike>(n: usize, g: &[C]) -> CharTableSummary {
    let mut ans = vec![];
    let glen = g.len();
    // add trivial repr
    {
        ans.push(CharEntry::new("triv".to_string(), 1.into(), &g));
    }

    if n >= 2 {
        // add alt repr
        let mut alt = CharEntry::new("alt".to_string(), 1.into(), &g);
        for i in 0..glen {
            let g = g[i].to_repr();
            if g.sgn() == -1 {
                alt.table[i] = (-1).into();
            }
        }
        ans.push(alt);
    }

    // add standard repr
    {
        let mut std = CharEntry::new("std".to_string(), 0.into(), &g);
        for i in 0..glen {
            let g = &g[i].to_repr();
            let fix = (0..n).filter(|&i| g[i] == i).count();
            std.table[i] = BigInt::from(fix) - 1;
        }
        ans.push(std);
    }

    let mut remaining = num_bigint::BigInt::from(glen);

    for r in &ans {
        // Assuming g[0] = e
        let dim = r.dim();
        remaining -= BigInt::from(&dim * &dim);
    }

    // try to find all representations, by creating tensor representations
    while remaining != 0.into() {
        if let Some(new_irrep) = find_new_entry(g, &ans) {
            // A new irrep was found.
            let dim = &new_irrep.table[0];
            remaining -= dim * dim;
            ans.push(new_irrep);
        } else {
            break;
        }
    }

    CharTableSummary {
        table: ans,
        remaining,
    }
}