use clap::Parser;
use rust_htslib::faidx::Reader;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};

/// Search for a specific sequence in the human pangenome
#[derive(Parser, Debug)]
struct Args {
    /// Path to the FASTA file
    #[arg(short, long)]
    fasta: String,

    /// Target sequence to search for
    #[arg(short, long)]
    sequence: String,

    /// Prefix for the sequence names to search within
    #[arg(short = 'p', long, default_value = "")]
    prefix: String,

    /// Maximum number of mismatches allowed (Hamming distance)
    #[arg(short, long, default_value_t = 6)]
    distance: usize,
}

fn load_fai(path: &str) -> HashMap<String, usize> {
    let fai_path = format!("{}.fai", path);
    let file = File::open(&fai_path).unwrap();
    let reader = BufReader::new(file);

    let mut sequences = HashMap::new();

    for line in reader.lines() {
        let line = line.unwrap();
        let parts: Vec<&str> = line.split_whitespace().collect();
        let name = parts[0].to_string();
        let length: usize = parts[1].parse().unwrap();
        sequences.insert(name, length);
    }

    sequences
}

fn hamming_distance(s1: &str, s2: &str) -> usize {
    s1.chars().zip(s2.chars()).filter(|&(c1, c2)| c1 != c2).count()
}

fn search_sequence(fasta: &str, target: &str, prefix: &str, max_mismatches: usize)  {
    let reader = Reader::from_path(fasta).unwrap();
    let n_seqs = reader.n_seqs();
    let seq_lengths = load_fai(fasta);
    // print a header line in tsv
    println!("seq_name\tstart\tend\tsequence\tmismatches");
    for i in 0..n_seqs {
        let seq_name = reader.seq_name(i as i32).unwrap();

        // Skip sequences that don't match the prefix
        if !prefix.is_empty() && !seq_name.starts_with(prefix) {
            continue;
        }
        let seq_length = seq_lengths.get(&seq_name).unwrap();
        let sequence_str = reader.fetch_seq_string(&seq_name, 0, *seq_length).unwrap();
        //println!("Searching sequence {} of length {}", seq_name, sequence_str.len());

        for (idx, window) in sequence_str.as_bytes().windows(target.len()).enumerate() {
            let window_str = std::str::from_utf8(window).unwrap();
            let distance = hamming_distance(window_str, target);
            if distance <= max_mismatches {
                //println!("Match found in sequence {}: {}", seq_name, window_str);
                println!("{}\t{}\t{}\t{}\t{}", seq_name, idx, idx+window.len(), window_str, distance);
            }
        }
    }
}

fn main() {
    let args = Args::parse();

    search_sequence(&args.fasta, &args.sequence, &args.prefix, args.distance);

}
