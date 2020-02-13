import os
import unittest

from jphmm_tools.subtyper import jphmm2bitmask, parse_bitmask, get_length, get_gappy_breakpoints
from jphmm_tools import save_bitmask, get_gap_mask
from hashlib import md5


DATA_DIR = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'data')
BREAKPOINTS = os.path.join(DATA_DIR, 'HIV1.breakpoints')
GAPPY_BREAKPOINTS = os.path.join(DATA_DIR, 'HIV1.gappy_breakpoints')
GAPMASK = os.path.join(DATA_DIR, 'HIV1.gapmask')
JPHMM_BITMASK = os.path.join(DATA_DIR, 'HIV1.bitmask')
ALN = os.path.join(DATA_DIR, 'HIV1.aln.fasta')
MSA = os.path.join(DATA_DIR, 'aln_to_msa.txt')
REC = os.path.join(DATA_DIR, 'recombination.txt')


class BitmaskTest(unittest.TestCase):

    def test_gappy_breakpoint_md5hash(self):
        out_bp = '{}.temp'.format(GAPPY_BREAKPOINTS)

        save_bitmask(get_gappy_breakpoints(breakpoint_file=BREAKPOINTS, aln_file=ALN)[0], out_bp)
        with open(GAPPY_BREAKPOINTS, 'r') as f:
            expected_hash = md5(f.read().encode()).hexdigest()
        with open(out_bp, 'r') as f:
            our_hash = md5(f.read().encode()).hexdigest()
        try:
            os.remove(out_bp)
        except:
            pass
        self.assertEqual(expected_hash, our_hash, msg='Hashsum of the file does not correspond to the expected one.')

    def test_gapmask(self):
        out_gm = '{}.temp'.format(GAPMASK)

        save_bitmask({name: {'-': gm} for name, gm in get_gap_mask(parse_bitmask(GAPPY_BREAKPOINTS)).items()}, out_gm)
        with open(GAPMASK, 'r') as f:
            expected_hash = md5(f.read().encode()).hexdigest()
        with open(out_gm, 'r') as f:
            our_hash = md5(f.read().encode()).hexdigest()
        try:
            os.remove(out_gm)
        except:
            pass
        self.assertEqual(expected_hash, our_hash, msg='Hashsum of the file does not correspond to the expected one.')

    def test_aln_len(self):
        crf2st2bm = parse_bitmask(GAPPY_BREAKPOINTS)
        n = get_length(crf2st2bm)
        self.assertEqual(35, n, msg='Got a wrong alignment length from gappy breakpoints')

    def test_jphmm_md5hash(self):
        crf2st2bm = parse_bitmask(GAPPY_BREAKPOINTS)
        n_gappy = get_length(crf2st2bm)
        id2insertion_length, id2n, id2st2bitmask = jphmm2bitmask(crf2st2bm, MSA, REC, n_gappy)
        out_bm = '{}.temp'.format(JPHMM_BITMASK)
        save_bitmask(id2st2bitmask, out_bm)
        with open(JPHMM_BITMASK, 'r') as f:
            expected_hash = md5(f.read().encode()).hexdigest()
        with open(out_bm, 'r') as f:
            our_hash = md5(f.read().encode()).hexdigest()
        try:
            os.remove(out_bm)
        except:
            pass
        self.assertEqual(expected_hash, our_hash, msg='Hashsum of the file does not correspond to the expected one.')

    def test_jphmm_lengths(self):
        crf2st2bm = parse_bitmask(GAPPY_BREAKPOINTS)
        n_gappy = get_length(crf2st2bm)
        id2insertion_length, id2n, id2st2bitmask = jphmm2bitmask(crf2st2bm, MSA, REC, n_gappy)
        for name, seq in [('nonCRF', 'AATTTTTTAAAAATTTTTTTTTTTTGG'), ('CRF2', 'GGGAAAAACCCCCCCCCCCCCCCAAAAA'), ('HXB2_or_CRF', 'AAAAAAAAAA')]:
            self.assertEqual(len(seq), id2n[name], msg='Sequence length miscalculated for {}.'.format(name))

    def test_jphmm_insertion_lengths(self):
        crf2st2bm = parse_bitmask(GAPPY_BREAKPOINTS)
        n_gappy = get_length(crf2st2bm)
        id2insertion_length, id2n, id2st2bitmask = jphmm2bitmask(crf2st2bm, MSA, REC, n_gappy)
        for name, n in [('nonCRF', 2), ('CRF2', 3), ('HXB2_or_CRF', 0)]:
            self.assertEqual(n, id2insertion_length[name], msg='Insertion length miscalculated for {}.'.format(name))