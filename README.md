# rust-schreier-sims [![Build Status](https://github.com/koba-e964/rust-schreier-sims/actions/workflows/rust.yml/badge.svg?branch=master)](https://github.com/koba-e964/rust-schreier-sims/actions/workflows/rust.yml?query=branch%3Amaster)
An implementation of Schreier-Sims algorithm in group theory.

This repository gives an implementation of [Schreier-Sims algorithm](https://en.wikipedia.org/wiki/Schreier%E2%80%93Sims_algorithm).
You can try this by executing main() function by `cargo run --release`:
```
fn main() {
    // Rubik's Cube group
    let n = 48;
    let mut gen = vec![];
    let p0 = vec![vec![0, 2, 4, 6], vec![1, 3, 5, 7],
                  vec![14, 16, 34, 28], vec![13, 23, 33, 27],
                  vec![12, 22, 32, 26]];
    let p2 = vec![vec![16, 18, 20, 22], vec![17, 19, 21, 23],
                  vec![12, 40, 36, 4], vec![11, 47, 35, 3],
                  vec![10, 46, 34, 2]];
    let p1 = vec![vec![8, 10, 12, 14], vec![9, 11, 13, 15],
                  vec![42, 18, 2, 26], vec![41, 17, 1, 25],
                  vec![40, 16, 0, 24]];
    let p3 = vec![vec![24, 26, 28, 30], vec![25, 27, 29, 31],
                  vec![8, 0, 32, 44], vec![15, 7, 39, 43],
                  vec![14, 6, 38, 42]];
    let p4 = vec![vec![32, 34, 36, 38], vec![33, 35, 37, 39],
                  vec![6, 22, 46, 30], vec![5, 21, 45, 29],
                  vec![4, 20, 44, 28]];
    let p5 = vec![vec![40, 42, 44, 46], vec![41, 43, 45, 47],
                  vec![10, 24, 38, 20], vec![9, 31, 37, 19],
                  vec![8, 30, 36, 18]];
    for p in &[p0, p1, p2, p3, p4, p5] {
        gen.push(get_cycle(n, p));
    }
    // Outputs order = 43252003274489856000
    eprintln!("order = {}", order(n, &gen));
}
```
In this example, `gen` becomes a set of generators of Rubik's Cube group.
The order of Rubik's Cube group is known to be 43252003274489856000 = 2<sup>11</sup> * 3<sup>7</sup> * 12! * 8! / 2, and this fact is confirmed by executing this program!
