# seqfu-amplify

A Python tool for performing in silico PCR analysis on FASTA sequences. 
This tool identifies potential amplicon regions in DNA sequences based on provided forward and reverse primers, 
based on [seqfu grep](https://telatin.github.io/seqfu2/tools/grep.html).

[**SeqFu**](https://telatin.github.io/seqfu2/) supports degenerate bases and mismatches. 

:warning: This is a preliminary implementation

## Requirements

- Python 3.6+
- SeqFu (must be installed and accessible in PATH)

## Installation

Clone this repository and ensure you have SeqFu installed:

```bash
# Install SeqFu if not already installed
conda install -c conda-forge -c bioconda python=3.9 "seqfu>=1.20"
```

## Usage

```bash
python amplifu.py [-h] [-f FWD_PRIMER] [-r REV_PRIMER] 
                        [-m MIN_LEN] [-x MAX_LEN] [-s] 
                        [--verbose] FASTA_FILES [FASTA_FILES ...]
```

### Arguments

- `-f`, `--fwd-primer`: Forward primer sequence (default: `CCTACGGGNGGCWGCAG`)
- `-r`, `--rev-primer`: Reverse primer sequence (default: `GGACTACHVGGGTATCTAATCC`)
- `-m`, `--min-len`: Minimum amplicon length (default: 100)
- `-x`, `--max-len`: Maximum amplicon length (default: 10000)
- `-s`, `--forward-strand`: Output the sequence in the forward strand
- `--verbose`: Print verbose output
- `FASTA_FILES`: One or more input FASTA files (can be gzipped)

### Output

The tool outputs a tab-separated file with the following columns:

- FileName: Name of the input FASTA file
- SeqName: Name of the sequence
- StartPosition: Start position of the amplicon
- EndPosition: End position of the amplicon
- FirstPrimer: Type of the first primer match
- SecondPrimer: Type of the second primer match
- AmpliconLength: Length of the amplicon
- AmpliconSequence: Sequence of the amplicon

## Example output

```bash
python amplifu.py --min-len 10 test/test.fa
```

| FileName | SeqName | StartPosition | EndPosition | FirstPrimer | SecondPrimer | AmpliconLength | AmpliconSequence |
|----------|----------|---------------|-------------|-------------|--------------|----------------|-----------------|
| test.fa | Ampli1_10 | 17 | 27 | forward_plus | reverse_minus | 10 | aaaaannnng |
| test.fa | Ampli1_10_rc | 22 | 32 | reverse_plus | forward_minus | 10 | AAAANNNNNT |
| test.fa | Ampli2_10 | 17 | 27 | forward_plus | reverse_minus | 10 | aaaaannnng |
| test.fa | Ampli2_10_rc | 22 | 32 | reverse_plus | forward_minus | 10 | AAAANNNNNT |
| test.fa | Ampli3_20_0 | 17 | 37 | forward_plus | reverse_minus | 20 | aaaaaNNNNNNNNNnnnnng |
| test.fa | Ampli3_20_0_rc | 22 | 42 | reverse_plus | forward_minus | 20 | AAAANNNNNNNNNNNNNNNT |
| test.fa | Ampli4_20_0 | 17 | 37 | forward_plus | reverse_minus | 20 | aaaaaNNNNNNNNNnnnnng |
| test.fa | Ampli4_20_0_rc | 22 | 42 | reverse_plus | forward_minus | 20 | AAAANNNNNNNNNNNNNNNT |
| test.fa | Ampli5_20_10 | 27 | 47 | forward_plus | reverse_minus | 20 | aaaaaNNNNNNNNNnnnnng |
| test.fa | Ampli5_20_10_rc | 22 | 42 | reverse_plus | forward_minus | 20 | AAAANNNNNNNNNNNNNNNT |
| test.fa | Ampli6_20_10 | 17 | 37 | forward_plus | reverse_minus | 20 | aaaaaNNNNNNNNNnnnnng |
| ...     |         |    |    |              |               |    |   |
| test.fa | Ampli_1_3 | 66 | 86 | forward_plus | reverse_minus | 20 | aaaaaNNNNNNNNNnnnnng |

 