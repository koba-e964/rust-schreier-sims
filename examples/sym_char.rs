use rust_schreier::sym_char::classes::get_all_concjugacy_classes;
use rust_schreier::sym_char::{get_char_table, CharTableSummary};

fn main() {
    let n = 10;
    let all = get_all_concjugacy_classes(n);
    let CharTableSummary {
        table: ans,
        remaining,
    } = get_char_table(n, &all);
    for ans in &ans {
        print!("{}:", ans.name);
        for j in 0..1 {
            print!("\t{}", ans.table[j]);
        }
        println!();
    }
    for i in 0..ans.len() {
        for j in 0..ans.len() {
            let prod = ans[i].inner_prod(&ans[j], &all);
            print!(" {}", prod);
        }
        println!();
    }
    println!("remaining = {}", remaining);
}
