import os
import unittest

from jphmm_tools.subtyper import get_subtypes

DATA_DIR = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'data')
FA = os.path.join(DATA_DIR, 'HIV1.fasta')
ALN = os.path.join(DATA_DIR, 'HIV1.fasta.aligned')
MSA = os.path.join(DATA_DIR, 'aln_to_msa.txt')
GAPPY_BREAKPOINTS = os.path.join(DATA_DIR, 'HIV1.gappy_breakpoints')
REC = os.path.join(DATA_DIR, 'recombination.txt')


class SubtypingTest(unittest.TestCase):

    def test_jphmm_subtypes(self):
        df = get_subtypes(MSA, REC, GAPPY_BREAKPOINTS)
        for name, st in [('CRF1', 'A,B'), ('nonCRF', 'A,B'), ('CRF2', 'A,C'), ('HXB2_or_CRF', 'B'), ('CRF', 'A,C')]:
            self.assertEqual(st, df.loc[name, 'subtype_jpHMM'], msg='Wrong subtype(s) for {}.'.format(name))

    def test_jphmm_compatible_subtypes(self):
        df = get_subtypes(MSA, REC, GAPPY_BREAKPOINTS)
        for name, st in [('CRF1', 'CRF1'), ('nonCRF', ''), ('CRF2', 'CRF2'), ('HXB2_or_CRF', 'B/CRF1'), ('CRF', 'CRF2')]:
            self.assertEqual(st, df.loc[name, 'compatible_subtypes'], msg='Wrong compatible subtype(s) for {}.'.format(name))