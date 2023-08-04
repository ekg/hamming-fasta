# hamming-fasta

Brute force sequence similarity search.

`hamming-fasta` searches for a specific sequence in a FASTA file, allowing for a specified number of mismatches.
It outputs any matches found along with their location in the genome and number of mismatches.
It can run in parallel and is designed to for "online" brute force sequence matching problems, such as where you want to find all occurrences of a particular string.
One application of this is CRISPR guide design.

## Usage

The program accepts the following arguments:

- `--fasta`: Path to the FASTA file to search in 
- `--sequence`: The sequence to search for
- `--prefix`: Only search sequences starting with this prefix (optional)  
- `--distance`: Maximum number of mismatches allowed (default: 6)
- `--parallelism`: Number of threads for parallel execution (default: number of CPU cores)

For example:

```
sequence_search --fasta genome.fa --sequence ACGT --distance 2
```

This will search genome.fa for the sequence ACGT, allowing up to 2 mismatches.

## Implementation

- The program first loads a lookup table of sequence names and lengths from the FASTA index file.

- It then iterates through the sequences in parallel, extracting each subsequence window the size of the target sequence. 

- The Hamming distance is calculated between the window and target sequence and matches below the mismatch threshold are printed.

- Access to standard output is synchronized via a mutex to avoid interleaved output.

- The rayon crate is used for simple parallelization across sequences.

## Output

The output is tab-separated with the following columns:

```
seq_name    start   end     sequence    mismatches
```

Where `seq_name` is the FASTA sequence name, `start` and `end` define the range of the match, `sequence` is the extracted sequence, and `mismatches` is the number of differences from the search sequence.
