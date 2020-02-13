import os
import unittest
from hashlib import md5

from jphmm_tools.aligner import align

DATA_DIR = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'data')
FA = os.path.join(DATA_DIR, 'HIV1.fasta')
ALN = os.path.join(DATA_DIR, 'HIV1.fasta.aligned')
MSA = os.path.join(DATA_DIR, 'aln_to_msa.txt')


class AlignmentTest(unittest.TestCase):

    def test_md5hash(self):
        out_aln = '{}.temp'.format(ALN)
        align(MSA, FA, out_aln, 35)
        with open(ALN, 'r') as f:
            expected_hash = md5(f.read().encode()).hexdigest()
        with open(out_aln, 'r') as f:
            our_hash = md5(f.read().encode()).hexdigest()
        try:
            os.remove(out_aln)
        except:
            pass
        self.assertEqual(expected_hash, our_hash, msg='Hashsum of the file does not correspond to the expected one.')