use clap::Parser;
use rust_htslib::faidx::Reader;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::io::Write;
use std::sync::{Arc, Mutex};
use rayon::prelude::*;

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

    /// Number of threads for parallel execution
    #[arg(short = 't', long, default_value = "0")]
    parallelism: usize,

    /// Flag for the presence of a Cas9 PAM sequence
    #[arg(short = 'c', long, default_value = "false")]
    cas9: bool,
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

fn get_reverse_complement(seq: &str) -> String {
    let mut rev_comp = String::new();
    for nucleotide in seq.chars().rev() {
        let complement = match nucleotide {
            'A' | 'a' => 'T',
            'T' | 't' => 'A',
            'G' | 'g' => 'C',
            'C' | 'c' => 'G',
            _ => 'N',
        };
        rev_comp.push(complement);
    }
    rev_comp
}

fn search_sequence(fasta: &str, target: &str, prefix: &str, max_mismatches: usize, cas9: bool) {
    let reader = Reader::from_path(fasta).unwrap();
    let n_seqs = reader.n_seqs();
    let seq_lengths = load_fai(fasta);

    let stdout_lock = Arc::new(Mutex::new(std::io::stdout()));

    let mut target = target.to_string();
    if cas9 {
        // add the PAM sequence to the target
        target.push_str("NGG");
        println!("seq_name\tstrand\tstart\tend\tsequence\tmismatches.nonpam");
    } else {
        // print a header line in tsv
        println!("seq_name\tstrand\tstart\tend\tsequence\tmismatches");
    }

    (0..n_seqs).into_par_iter().for_each(|i| {
        let reader = Reader::from_path(fasta).unwrap(); // Re-create the reader for thread safety
        let seq_name = reader.seq_name(i as i32).unwrap();
        if !prefix.is_empty() && !seq_name.starts_with(prefix) {
            return;
        }
        let seq_length = seq_lengths.get(&seq_name).unwrap();
        let sequence_str = reader.fetch_seq_string(&seq_name, 0, *seq_length).unwrap();
        let rev_sequence_str = get_reverse_complement(&sequence_str);

        for sequence in [(sequence_str, "+"), (rev_sequence_str, "-")].iter() {
            for (idx, window) in sequence.0.as_bytes().windows(target.len()).enumerate() {
                let window_str = std::str::from_utf8(window).unwrap();
                if cas9 && !window_str.ends_with("GG") {
                    continue;
                }
                let mut distance = hamming_distance(window_str, &target);
                if cas9 {
                    // because we have N in the query sequence, we need to subtract 1 from the distance
                    distance -= 1;
                }
                if distance <= max_mismatches {
                    let mut stdout = stdout_lock.lock().unwrap();

                    let start = if sequence.1 == "-" { *seq_length - (idx + window.len()) } else { idx };
                    let end = if sequence.1 == "-" { *seq_length - idx } else { idx + window.len() };

                    writeln!(stdout, "{}\t{}\t{}\t{}\t{}\t{}", seq_name, sequence.1, start, end, window_str, distance).unwrap();
                }
            }
        }
    });
}

fn main() {
    let args = Args::parse();
    // If parallelism is set to 0, use the default (number of available CPU cores)
    if args.parallelism > 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(args.parallelism)
            .build_global()
            .unwrap();
    }
    search_sequence(&args.fasta, &args.sequence, &args.prefix, args.distance, args.cas9);
}
