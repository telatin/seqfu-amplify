#!/usr/bin/env python3
"""
Test suite for in silico PCR analysis
"""
import unittest
import tempfile
import os
from pathlib import Path
import subprocess

class TestInSilicoPCR(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        """Create test data files"""
        cls.test_dir = tempfile.mkdtemp()
        cls.script_path = Path("amplifu.py")
        
        # Create test FASTA file
        cls.test_fasta = os.path.join(cls.test_dir, "test.fa")
        with open(cls.test_fasta, "w") as f:
            f.write(""">Ampli1_10
NNNNNNNNNNNCAGATANNNNNNNNNNNNNGGTTTGGNNNNNNNNNN
>Ampli1_10_rc
NNNNNNNNNNCCAAAACCNNNNNNNNNNNNNNTATCTGNNNNNNNNNN
>Ampli2_10
NNNNNNNNNNNCAGATANNNNNNNNNNNNNGGTTTTGGnnnnnnnnnn
>Ampli2_10_rc
nnnnnnnnnnCCAAAAAACNNNNNNNNNNNNNTATCTGnnnnnnnnnn
>Ampli3_20_0
NNNNNNNNNNCAGATANNNNNNNNNNNNNNGGTTTTGG
>Ampli3_20_0_rc
CCAAAACCNNNNNNNNNNNNNNTATCTGNNNNNNNNNN
>Ampli_1_3
NNNNNNNNNNCAGATANNNNNNNNNNNNNGGTTTGGNNNGGATTAGATACCCTGGTAGTCCCCCTACGGGTGGCAGCAGNNNNCAGATANNNNNNNNNNNNNGGTTTGGNNNN
""")
        
    def test_basic_amplification(self):
        """Test basic amplification with default parameters"""
        cmd = [
            "python3",
            str(self.script_path),
            "-f", "CAGATA",
            "-r", "CCAAACC",
            "-m", "10",
            self.test_fasta
        ]
        result = subprocess.run(cmd, capture_output=True, text=True)
        self.assertEqual(result.returncode, 0)
        
        # Parse output and check results
        lines = result.stdout.strip().split("\n")
        headers = lines[0].split("\t")
        data = [line.split("\t") for line in lines[1:]]
        
        # Check that we found the expected number of amplicons
        self.assertGreater(len(data), 0, "No amplicons found")
        
        # Verify structure of output
        self.assertEqual(len(headers), 8, "Incorrect number of columns in output")
        expected_headers = [
            "FileName", "SeqName", "StartPosition", "EndPosition",
            "FirstPrimer", "SecondPrimer", "AmpliconLength", "AmpliconSequence"
        ]
        self.assertEqual(headers, expected_headers)

    def test_size_constraints(self):
        """Test amplicon size constraints"""
        # Test with very small max length
        cmd = [
            "python3",
            str(self.script_path),
            "-f", "CAGATA",
            "-r", "CCAAACC",
            "-x", "15",  # Max length 15
            self.test_fasta
        ]
        result = subprocess.run(cmd, capture_output=True, text=True)
        self.assertEqual(result.returncode, 0)
        
        # All amplicons should be ≤ 15bp
        lines = result.stdout.strip().split("\n")[1:]  # Skip header
        for line in lines:
            length = int(line.split("\t")[6])  # AmpliconLength column
            self.assertLessEqual(length, 15)

    def test_strand_orientation(self):
        """Test forward strand output option"""
        cmd = [
            "python3",
            str(self.script_path),
            "-f", "CAGATA",
            "-r", "CCAAACC",
            "-s",  # Forward strand option
            self.test_fasta
        ]
        result = subprocess.run(cmd, capture_output=True, text=True)
        self.assertEqual(result.returncode, 0)
        
        # Check that reverse complement sequences are properly oriented
        lines = result.stdout.strip().split("\n")[1:]  # Skip header
        for line in lines:
            fields = line.split("\t")
            if "reverse" in fields[4]:  # If first primer is on reverse strand
                # Sequence should be in forward orientation
                self.assertFalse(fields[7].startswith("CCAAACC"))

    def test_invalid_input(self):
        """Test handling of invalid inputs"""
        # Test with non-existent file
        cmd = [
            "python3",
            str(self.script_path),
            "nonexistent.fa"
        ]
        result = subprocess.run(cmd, capture_output=True, text=True)
        self.assertNotEqual(result.returncode, 0)

    @classmethod
    def tearDownClass(cls):
        """Clean up test files"""
        try:
            os.remove(cls.test_fasta)
            os.rmdir(cls.test_dir)
        except:
            pass

if __name__ == "__main__":
    unittest.main(verbosity=2)
