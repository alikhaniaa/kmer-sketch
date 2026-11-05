//! kmer-sketch: canonical k-mer counting, FracMinHash sketching, sketch
//! comparison, read simulation and single-linkage clustering of cells.
//!
//! Suggested reading order: `kmer` -> `hash` -> `count` -> `sketch` ->
//! `compare` -> `simulate` -> `cluster`, then `main.rs` for the CLI.

pub mod cluster;
pub mod compare;
pub mod count;
pub mod hash;
pub mod kmer;
pub mod seqio;
pub mod simulate;
pub mod sketch;
