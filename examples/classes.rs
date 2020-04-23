use rust_schreier::sym_char::classes::{get_all_concjugacy_classes, ClassLike};

fn main() {
    let n = 6;
    let all = get_all_concjugacy_classes(n);
    for v in &all {
        let repr = v.to_repr();
        println!("{:?}, {}", repr, v.get_size());
    }
}
