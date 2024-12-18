#!/usr/bin/env python3
"""
>ampli_nopad
CAGATANNNNNNNNNNNNNNGGTTTTGG
>ampli_pad10
nnnnnnnnnnCAGATANNNNNNNNNNNNNNGGTTTTGGnnnnnnnnnn
>revampli_nopad
CCAAAACC18+14NNNNNNNNNNNNNNTATCTG
>revampli_pad10
NNNNNNNNNNCCAAAACCNNNNNNNNNNNNNNTATCTGNNNNNNNNNN

"""
import argparse
import re
import subprocess
import sys
from typing import List, Tuple, Dict
from pathlib import Path
import gzip

def fasta_generator(filename):
    """
    Parse a FASTA (gzipped or not) file and yield header and sequence (id, comments, seq)
    """
    with gzip.open(filename, 'rt') if filename.endswith('.gz') else open(filename) as f:
        sequence = []
        seqid, comments = None, None
        for line in f:
            if line.startswith('>'):
                if sequence:
                    seqid, *comments = header[1:].split()
                    comments = ' '.join(comments)
                    yield seqid, comments, ''.join(sequence)
                header = line.strip()
                seqid, *comments = header[1:].split()
                comments = ' '.join(comments)
                sequence = []
            else:
                sequence.append(line.strip())
        if sequence:
            yield seqid, comments, ''.join(sequence)

def get_slice(sequence, start, end):
    return sequence[start:end]

def check_seqfu():
    cmd = ['seqfu', '--version']
    version = None
    try:
        execution = subprocess.run(cmd, capture_output=True, text=True, check=True)
        version = execution.stdout.strip()
    except subprocess.CalledProcessError as e:
        print(f"Error running seqfu: {e}")
    return version

def find_amplicons(primerFor_matchesPlus, primerFor_matchesMinus, primerRev_matchesPlus, primerRev_matchesMinus, min_len, max_len):
    """
    Find amplicons from primer matches.
    Get the start positions of FORWARD primer (in both plus and minus strands) and REVERSE primer (in both plus and minus strands)

    Valid amplicons: 
    1. [FWDPRIM_plus -> REVPRIM_minus]: Forward primer on plus strand to Reverse primer on minus strand
    2. [REVPRIM_plus -> FWDPRIM_minus]: Reverse primer on plus strand to Forward primer on minus strand
    
    Args:
        primerFor_matchesPlus: List of positions where forward primer matches on plus strand
        primerFor_matchesMinus: List of positions where forward primer matches on minus strand
        primerRev_matchesPlus: List of positions where reverse primer matches on plus strand
        primerRev_matchesMinus: List of positions where reverse primer matches on minus strand
        min_len: Minimum allowed amplicon length
        max_len: Maximum allowed amplicon length
        
    Returns:
        List of tuples (start, end, first_primer, second_primer)
    """
    amplicons = []
    
    # Case 1: Forward primer (plus strand) to Reverse primer (minus strand)
    # Direction: 5' -> 3' on plus strand
    for start in primerFor_matchesPlus:
        for end in primerRev_matchesMinus:
            length = end - start
            if min_len <= length <= max_len:
                amplicons.append((start, end, "forward_plus", "reverse_minus"))
    
    # Case 2: Reverse primer (plus strand) to Forward primer (minus strand)
    # Direction: 5' -> 3' on minus strand
    for start in primerRev_matchesPlus:
        for end in primerFor_matchesMinus:
            length = end - start
            if min_len <= length <= max_len:
                amplicons.append((start, end, "reverse_plus", "forward_minus"))
    
    return sorted(amplicons, key=lambda x: (x[0], x[1]))

def run_seqfu_grep(oligo: str, fasta_file: str, verbose=False) -> str:
    """Run seqfu grep command and return output."""
    cmd = ['seqfu', 'grep', '-A', '-o', oligo, fasta_file]
    if verbose:
        print(f"# Running seqfu grep for {fasta_file} with command: {' '.join(cmd)}", file=sys.stderr)
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        lines = [line for line in result.stdout.splitlines() if line.startswith('>')]
    
        return lines
    except subprocess.CalledProcessError as e:
        print(f"Error running seqfu grep for {fasta_file}: {e}")
        return ""

def reverse_complement(seq: str) -> str:
    """Return reverse complement of sequence."""
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N', 'a': 't', 'c': 'g', 'g': 'c', 't': 'a', 'n': 'n'}
    # Check if there are invalid characters
    if not re.match('^[ACGTNacgtn]+$', seq):
        raise ValueError(f"Invalid sequence: {seq}")
    return ''.join(complement[base] for base in reversed(seq))

def parse_seqfu_output(output_lines: str) -> Dict[str, str]:
    """Parse seqfu output and return dictionary of sequence names and match lines."""
    matches = {}
    current_header = None
    
    for line in output_lines:
        current_header = line.split()[0][1:]  # Remove '>' and get first word
        # line contains "for-matches=:rev-matches=2012577,2590732", return this part
        for_parts = line.split('for-matches=')[1].split(':rev-matches=')[0]
        rev_parts = line.split(':rev-matches=')[1]

        # make int every item for for_parts and rev_parts
        for_parts = [int(x) for x in for_parts.split(',') if x]
        rev_parts = [int(x) for x in rev_parts.split(',') if x]
        matches[current_header] = (for_parts, rev_parts)
    
    return matches

def main():
    parser = argparse.ArgumentParser(description='In silico PCR analysis')
    parser.add_argument('-f', '--fwd-primer', default="CCTACGGGNGGCWGCAG", help='Forward primer sequence (default: %(default)s)')
    parser.add_argument('-r', '--rev-primer', default="GGACTACHVGGGTATCTAATCC", help='Reverse primer sequence (default: %(default)s)')
    parser.add_argument('-m', '--min-len', type=int, default=100,
                        help='Minimum amplicon length  (default: %(default)s)')
    parser.add_argument('-x', '--max-len', type=int, default=10000,
                        help='Maximum amplicon length  (default: %(default)s)')
    parser.add_argument('-s', '--forward-strand',  action='store_true', help='Output the sequence in the forward strand')
    parser.add_argument('--verbose', action='store_true', help='Print verbose output')
    parser.add_argument('fasta_files', nargs='+', help='Input FASTA files')
    
    args = parser.parse_args()
    
    verbose = args.verbose
    seqfu_version = check_seqfu()

    fwd_len = len(args.fwd_primer)
    rev_len = len(args.rev_primer)

    if verbose:
        print(f"# Using seqfu version: {seqfu_version}", file=sys.stderr)
    # Print header
    print("FileName\tSeqName\tStartPosition\tEndPosition\tFirstPrimer\tSecondPrimer\tAmpliconLength\tAmpliconSequence")
    
    for fasta_file in args.fasta_files:
        fasta_path = Path(fasta_file)
        
        sequence = {}
        for seqid, comments, seqstr in fasta_generator(fasta_file):
            sequence[seqid] = seqstr

        if verbose:
            print(f"# Processing {fasta_file} with {len(sequence)} sequences", file=sys.stderr)
        # Run seqfu grep for both primers
        fwd_output = run_seqfu_grep(args.fwd_primer, fasta_file, verbose)
        rev_output = run_seqfu_grep(args.rev_primer, fasta_file, verbose)
        
        # Parse outputs
        fwd_matches = parse_seqfu_output(fwd_output)
        rev_matches = parse_seqfu_output(rev_output)
 
        # Find sequences that match both primers
        common_seqs = set(fwd_matches.keys()) & set(rev_matches.keys())
        if verbose:
            print(f"# Found {len(common_seqs)} sequences with both primers", file=sys.stderr)
        # Sort common_seqs
        common_seqs = sorted(common_seqs)
        for seq_name in common_seqs:
            # Parse matches from headers
            primerFor_matchesPlus, primerFor_matchesMinus = fwd_matches[seq_name]
            primerRev_matchesPlus, primerRev_matchesMinus = rev_matches[seq_name]   
            
            # Find valid amplicons
            amplicons = find_amplicons(
                primerFor_matchesPlus, primerFor_matchesMinus, 
                primerRev_matchesPlus, primerRev_matchesMinus,
                args.min_len, args.max_len)
            if verbose:
                print(f"# Found {len(amplicons)} amplicons for {seq_name}", file=sys.stderr)
             
            # Output results
            for start, end, first_primer, second_primer in amplicons:
                if first_primer == "forward_plus" and second_primer == "reverse_minus":
                    start += fwd_len    
                elif first_primer == "reverse_plus" and second_primer == "forward_minus":
                    start += rev_len
                seqslice = get_slice(sequence[seq_name], start, end)
                if args.forward_strand:
                    if first_primer == "reverse_plus" and second_primer == "forward_minus":
                        seqslice = reverse_complement(seqslice)
                print(f"{fasta_path.name}\t{seq_name}\t{start}\t{end}\t{first_primer}\t{second_primer}\t{end-start}\t{seqslice}")

if __name__ == '__main__':
    main()