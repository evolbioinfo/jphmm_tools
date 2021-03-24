import logging

from Bio import SeqIO
from Bio.Alphabet import generic_dna

from jphmm_tools import VERSION


def split(in_fa, out_fa_pattern=None, chunk_size=1, n=0, exclude_ids=None):
    """
    Splits an input fasta file into several files.

    :param in_fa: path to the input fasta file
    :type in_fa: str
    :param out_fa_pattern: path pattern for output fasta files
        (should contain {} that will be replaced with the 0-based file number,
        e.g. "out.{}.fa" that will produce "out.0.fa", "out.1.fa", etc.).
        If not specified the in_fa value will be suffixed with 0-based numbers
    :type out_fa_pattern: str
    :param chunk_size: number of sequences to be placed into each split file
    :type chunk_size: int
    :param n: total number of output split files. If specified, the chunk_size will be adjusted accordingly
    :type n: int
    :param exclude_ids: names of the sequences in the input fasta that should be ignored and not put into the split files
    :type exclude_ids: iterable(str)
    """
    if not out_fa_pattern:
        out_fa_pattern = in_fa + '.{}'

    exclude_ids = {_.strip('\r').strip() for _ in exclude_ids} if exclude_ids else []
    if exclude_ids:
        logging.info('Excluding the sequences {} from the split files...'.format(exclude_ids))

    sequences = [rec for rec in SeqIO.parse(in_fa, 'fasta', alphabet=generic_dna) if rec.id not in exclude_ids]
    if n:
        chunk_size = len(sequences) // n
    elif chunk_size:
        if len(sequences) % chunk_size:
            n = len(sequences) // (chunk_size - 1)
        else:
            n = len(sequences) // chunk_size
    else:
        n = len(sequences)
        chunk_size = 1
    remainder = len(sequences) % n

    logging.info('Splitting {} sequences into {} files...'.format(len(sequences), n))
    i = 0
    start = 0
    while start < len(sequences):
        stop = start + chunk_size
        if i < remainder:
            stop += 1
        SeqIO.write(sequences[start: min(len(sequences), stop)], out_fa_pattern.format(i), 'fasta')
        start = stop
        i += 1


def main():
    """
    Entry point, calling :py:func:`jphmm_tools.aln_splitter.split` with command-line arguments.

    :return: void
    """
    import argparse

    parser = argparse.ArgumentParser(description="Splits an input fasta file into several files.", prog='jphmm_split')

    parser.add_argument('--version', action='version', version='%(prog)s {version}'.format(version=VERSION))

    parser.add_argument('--in_fa', required=True, type=str, help='Path to the input fasta file.')
    parser.add_argument('--out_fa_pattern', required=False, type=str,
                        help='Path pattern for output fasta files '
                             '(should contain {} that will be replaced with the 0-based file number, '
                             'e.g. "out.{}.fa" that will produce "out.0.fa", "out.1.fa", etc.). '
                             'If not specified the --in_fa value will be suffixed with 0-based numbers.')
    parser.add_argument('--chunk_size', required=False, type=int, default=1,
                        help="Number of sequences to be placed into each split file.")
    parser.add_argument('--n', required=False, type=int, default=None,
                        help="Total number of output split files. "
                             "If specified, the --chunk_size will be adjusted accordingly.")
    parser.add_argument('--exclude_ids', nargs='*', type=str,
                        help="Names of the sequences in the input fasta that should be ignored "
                             "and not put into the split files.")
    params = parser.parse_args()

    logging.basicConfig(level=logging.INFO)
    split(params.in_fa, params.out_fa_pattern, chunk_size=params.chunk_size, n=params.n,
          exclude_ids=params.exclude_ids)


if '__main__' == __name__:
    main()
