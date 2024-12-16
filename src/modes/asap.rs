//! # The rusty asap
//!
//! The main purpose of the ASAP algorithm is to find the hclust
//! cutoff value that is optimal for seperating amplicons into species.
//!
//! The asap algorithm has the following steps:
//! 1. Ascending hierarchical clustering.
//! 2. Computing P values
//! 3. Recursive splits
//! 4. Relative barcode gap width.

use crate::core::analysis::distance::DistanceAnalysisBuilder;
use crate::core::io;
use crate::core::utils::pairs::{Pair, Pairwise};
use crate::core::utils::{fasta_distance_jukes_cantor_number, remove_empty};
use bio::io::fasta::Record;
use clap::Args;
use rayon::iter::IntoParallelIterator;
use rayon::iter::ParallelIterator;
use rayon::prelude::IntoParallelRefIterator;
use std::error::Error;
use std::fmt::format;

fn find_length_differences(mut input: Vec<Record>) -> Vec<Pair<Record, f64>> {
    fn len_diff(
        x: Record,
        y: Record,
    ) -> Pair<Record, f64> {
        let diff = x.seq().len().abs_diff(y.seq().len()) as f64;
        Pair::new(x, y, diff)
    }
    let len = input.pairwise_map_condensed_upper(len_diff);
    len.into_par_iter().filter(|x| x.x != 0_f64).collect()
}

// verify
fn verify_records(
    records: Vec<Record>
) -> Result<Vec<Record>, Box<dyn Error>> {
    if records.len() < 2 {
        Err(Box::from("There are only 2 records in the file! More records are needed to delineate species."))
    } else {
        let lengths = records
            .iter()
            .cloned()
            .map(|x: Record| x.seq().len())
            .collect::<Vec<usize>>();
        let first = lengths[0];
        let differences = find_length_differences(records.clone());
        let max_diff =
            differences.iter().max_by(|x, y| x.partial_cmp(y).unwrap());
        match max_diff {
            Some(diff) => {
                let amount = records.len() - differences.len();
                let error = format!("There are {} differences in length in the file ({} records total). For instance, {} and {} differ",
                                    amount,
                                    records.len(),
                                    diff.a.id(),
                                    diff.b.id()
                );
                Err(Box::from(error))
            }
            None => Ok(records),
        }
    }
}

/// Gives back an Array of Fasta records valid for analysis.
///
/// Empty records are removed. If there are only two records or fewer. The function raises an error.
fn take_input(input: &str) -> Result<Vec<Record>, Box<dyn Error>> {
    let fasta = io::read_fasta(&input)?;
    let records = remove_empty(fasta);
    let verified = verify_records(records)?;
    Ok(verified)
}

/// Compute pairwise distances
fn compute_distances(mut records: Vec<Record>) -> Vec<Pair<Record, f64>> {
    let mut fasta_distances = records
        .pairwise_map_condensed_upper(fasta_distance_jukes_cantor_number);
    fasta_distances.sort_by(|a, b| a.partial_cmp(b).unwrap());
    fasta_distances
}

/// The entry point of the asap algorithm
pub(crate) fn asap(args: &str) -> Result<(), Box<dyn Error>> {
    let recs = take_input("resources/test/data/asv-listerria-taxon-Bacillales-Order.fasta.final_tree.fa")?;
    let recs = remove_empty(recs);
    let distanceAnalysis = DistanceAnalysisBuilder::new()
        .data(recs)
        .f(fasta_distance_jukes_cantor_number)
        .build()
        .unwrap();
    let mat = distanceAnalysis.run();
    let amount = mat.len();
    let top = mat.max();
    let distance = top.x;
    let msg = format!("{distance:?}");
    println!("Most different pair\n{:?}", top);
    println!("{} distances calculated", amount);
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::io::read_fasta;

    /// This test may give an IO error if the test file is not found.
    #[test]
    fn test_input() -> Result<(), Box<dyn Error>> {
        // File does not exist
        let recs =
            read_fasta("resources/test/data/sample.fasta-2line.un-eq-len")?;
        let recs = remove_empty(recs);
        let out = verify_records(recs);
        if let Err(err) = out {
            assert_eq!(err.to_string(), "There are 1 differences in length in the file (471 records total). For instance, 12a.1 and 293.1 differ");
        }

        // Correct file
        let recs = read_fasta("resources/test/data/asv-listerria-taxon-Bacillales-Order.fasta.final_tree.fa")?;
        let recs = remove_empty(recs);
        let out = verify_records(recs);
        match out {
            Err(err) => {
                panic!("This file should be correct!")
            }
            Ok(x) => assert_eq!(x.len(), 17),
        }

        // Correct file, too few sequences
        let recs = read_fasta("resources/test/data/too_few_records.fa")?;
        let recs = remove_empty(recs);
        let out = verify_records(recs);
        if let Err(err) = out {
            assert_eq!(err.to_string(), "There are only 2 records in the file! More records are needed to delineate species.");
        }

        Ok(())
    }

    #[test]
    fn run_asap() -> Result<(), Box<dyn Error>> {
        asap(
            "resources/test/data/asv-listerria-taxon-Bacillales-Order.fasta",
        )?;
        Ok(())
    }
}
