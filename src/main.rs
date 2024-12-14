//! The main module for ASAP

use clap::{arg, Parser};
mod lib;
mod modes;
use modes::asap;

/// rust clone of ASAP.
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Input data. For fasta files, a chosen distance metric is calculated.
    #[arg(name = "Input data")]
    input: String,
}

// You should have a yaml document for the parameters of the code.

fn main() {
    let args = Args::parse();
    let out = asap::asap(args);
    match out {
        Ok(_output) => {}
        Err(error) => {
            panic!("{}", error)
        }
    }
}
