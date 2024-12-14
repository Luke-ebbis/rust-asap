//! # Input and output of files
use bio::io::fasta;
use bio::io::fasta::FastaRead;
use fasta::Record;
use std::fs::File;
use std::io::Error;
use std::io::{BufReader, ErrorKind};
pub fn read_fasta(input: &str) -> Result<Vec<Record>, Error> {
    let f = File::open(input)?;
    let reader_file = BufReader::new(f);
    let mut reader = fasta::Reader::new(reader_file);
    let mut record = Record::new();
    let mut records: Vec<Record> = Vec::new();
    reader.read(&mut record).expect("Failed to parse record");
    while !record.is_empty() {
        let check = record.check();
        if check.is_err() {
            Err(Error::new(
                ErrorKind::InvalidData,
                "The FASTA has failed to parse.",
            ))?;
        }
        reader.read(&mut record).expect("Failed to parse record");
        records.push(record.clone());
    }
    Ok(records)
}
