//! Reading FASTA/FASTQ files (plain or gzip-compressed) with needletail, and
//! small helpers for file names.

use anyhow::{bail, Context, Result};
use needletail::errors::ParseErrorKind;
use std::path::Path;

/// Call `handle` once for every sequence in a FASTA or FASTQ file.
/// The format and gzip compression are detected automatically.
/// An empty file (also an empty `.gz`) counts as a file with zero records.
pub fn for_each_sequence(path: &Path, mut handle: impl FnMut(&[u8])) -> Result<()> {
    if !path.is_file() {
        bail!("{} is not a file (or does not exist)", path.display());
    }
    let mut reader = match needletail::parse_fastx_file(path) {
        Ok(reader) => reader,
        Err(e) if matches!(e.kind, ParseErrorKind::EmptyFile) => return Ok(()),
        Err(e) => {
            return Err(e).with_context(|| format!("cannot open sequence file {}", path.display()))
        }
    };
    while let Some(record) = reader.next() {
        let record =
            record.with_context(|| format!("bad record in sequence file {}", path.display()))?;
        handle(&record.seq());
    }
    Ok(())
}

/// Read a sequence file in batches of roughly `batch_bases` bases, and call
/// `handle` with each batch. Batching lets the caller process many short
/// reads in parallel without loading the whole file into memory.
pub fn for_each_batch(
    path: &Path,
    batch_bases: usize,
    mut handle: impl FnMut(&[Vec<u8>]),
) -> Result<()> {
    let mut batch: Vec<Vec<u8>> = Vec::new();
    let mut bases_in_batch = 0;
    for_each_sequence(path, |seq| {
        bases_in_batch += seq.len();
        batch.push(seq.to_vec());
        if bases_in_batch >= batch_bases {
            handle(&batch);
            batch.clear();
            bases_in_batch = 0;
        }
    })?;
    if !batch.is_empty() {
        handle(&batch);
    }
    Ok(())
}

/// Read every sequence of a file into memory (used for reference genomes).
pub fn read_all_sequences(path: &Path) -> Result<Vec<Vec<u8>>> {
    let mut sequences = Vec::new();
    for_each_sequence(path, |seq| sequences.push(seq.to_vec()))?;
    Ok(sequences)
}

/// File extensions that are removed to turn a file name into a sample name.
const SEQUENCE_EXTENSIONS: [&str; 7] =
    [".fastq", ".fq", ".fasta", ".fa", ".fna", ".sketch", ".txt"];

/// Turn a path like `cells/ecoli__cell001.fastq.gz` into a sample name
/// like `ecoli__cell001`: drop the directory, a `.gz` suffix, and one
/// known sequence extension.
pub fn sample_name(path: &Path) -> String {
    let file_name = path
        .file_name()
        .map(|name| name.to_string_lossy().into_owned())
        .unwrap_or_default();
    let mut name = file_name.strip_suffix(".gz").unwrap_or(&file_name);
    for extension in SEQUENCE_EXTENSIONS {
        if let Some(stripped) = name.strip_suffix(extension) {
            name = stripped;
            break;
        }
    }
    name.to_string()
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    #[test]
    fn sample_name_strips_directory_and_extensions() {
        let cases = [
            ("cells/ecoli__cell001.fastq.gz", "ecoli__cell001"),
            ("genomes/ecoli.fna", "ecoli"),
            ("x.fa", "x"),
            ("sketches/a.sketch", "a"),
            ("noext", "noext"),
        ];
        for (path, expected) in cases {
            assert_eq!(sample_name(&PathBuf::from(path)), expected);
        }
    }

    fn count_records(path: &Path) -> Result<usize> {
        let mut n = 0;
        for_each_sequence(path, |_| n += 1)?;
        Ok(n)
    }

    #[test]
    fn empty_files_have_zero_records() {
        use flate2::write::GzEncoder;
        use flate2::Compression;

        let dir = tempfile::tempdir().unwrap();
        let plain = dir.path().join("empty.fastq");
        std::fs::write(&plain, b"").unwrap();
        assert_eq!(count_records(&plain).unwrap(), 0);

        // A valid gzip stream with no content inside.
        let gz = dir.path().join("empty.fastq.gz");
        let encoder = GzEncoder::new(std::fs::File::create(&gz).unwrap(), Compression::default());
        encoder.finish().unwrap();
        assert_eq!(count_records(&gz).unwrap(), 0);
    }

    #[test]
    fn directory_and_garbage_give_errors() {
        let dir = tempfile::tempdir().unwrap();
        let err = count_records(dir.path()).unwrap_err();
        assert!(err.to_string().contains("is not a file"));

        let bad = dir.path().join("bad.fa");
        std::fs::write(&bad, b"this is not fasta\n").unwrap();
        assert!(count_records(&bad).is_err());
    }
}
