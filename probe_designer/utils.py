from Bio.Blast.Applications import NcbiblastnCommandline
from io import StringIO
import logging
import os
import pandas as pd


logging.basicConfig(format='%(levelname)s\t%(asctime)s\t%(name)s\t%(message)s',
                    level=logging.INFO)
logger = logging.getLogger('SeqUtils')


class FastaReader:
    def __init__(self, fasta, fasta_idx=None):
        """Reference fasta reader

        Parameters
        ----------
        fasta : path
            Path to reference fasta
        fasta_idx : path, optional
            Path to reference fasta index

        Raises
        ------
        FileNotFoundError
            Fasta index is required, raises error if not found
        """
        if not fasta_idx:
            fasta_idx = fasta + '.fai'

        if not os.path.exists(fasta_idx):
            raise FileNotFoundError('Fasta file must be indexed (.fasta.fai)')

        self.__fasta_idx_df = pd.read_table(
            fasta_idx, header=None, dtype={'contig': str},
            names=['contig', 'length', 'start',
                   'bases_per_line', 'bytes_per_line'])
        self.__fasta_idx_df.set_index('contig', inplace=True)
        self.__fasta_fh = open(fasta, 'r')
        self.fasta_path = fasta

    def get(self, contig, seq_start, seq_end):
        """Get sequence from fasta

        Parameters
        ----------
        contig : str
            Chromosome
        seq_start : int
            Start position [0-indexed]
        seq_end : int
            End position [1-indexed]

        Returns
        -------
        str
            Returns sequence string
        """
        contig = str(contig)
        if seq_start < 0:
            raise ValueError('Sequence start must be greater than 0')
        bases_per_line = self.__fasta_idx_df.loc[contig, 'bases_per_line']
        bytes_per_line = self.__fasta_idx_df.loc[contig, 'bytes_per_line']
        seq_length = seq_end - seq_start
        contig_start = self.__fasta_idx_df.loc[contig, 'start']
        num_newlines = seq_start // bases_per_line

        adj_seq_start = contig_start + seq_start + num_newlines
        max_seq_newlines = (seq_length // bytes_per_line) + 2
        max_seq_length = seq_length + max_seq_newlines

        self.__fasta_fh.seek(adj_seq_start)
        raw_seq = self.__fasta_fh.read(max_seq_length)
        raw_seq = raw_seq.replace('\n', '')
        seq = raw_seq[:seq_length]

        return seq

    def close(self):
        self.__fasta_fh.close()


def run_blast(sequence, reference, qcov_hsp_perc, perc_identity,
              threads=1, **kwargs):
    """ Run blast on a sequence and return results """
    if isinstance(sequence, str):
        fasta = f'>id\n{sequence}\n'
        length = len(sequence)
    elif isinstance(sequence, dict):
        fasta = '\n'.join(f'>{k}\n{v}' for k, v in sequence.items())
        length = len(next(iter(sequence.values())))

    blast_dtypes = {'qseqid': str, 'sseqid': str, 'pident': float,
                    'length': int, 'mismatch': int, 'gapopen': int,
                    'qstart': int, 'qend': int, 'sstart': int,
                    'send': int, 'evalue': str, 'bitscore': float,
                    'qcovhsp': float, 'qseq': str, 'sseq': str,
                    'sstrand': str}

    blast_cols = blast_dtypes.keys()
    blast_fmt = '"6 {}"'.format(' '.join(blast_cols))
    task = "blastn-short" if length <= 40 else "blastn"

    blastn_cline = NcbiblastnCommandline(db=reference,
                                         outfmt=blast_fmt,
                                         qcov_hsp_perc=qcov_hsp_perc,
                                         perc_identity=perc_identity,
                                         dust='no',
                                         soft_masking=False,
                                         task=task,
                                         num_threads=threads,
                                         **kwargs)
    out, err = blastn_cline(stdin=fasta)
    blast_df = pd.read_table(StringIO(out), names=blast_cols, dtype=blast_dtypes)
    return blast_df


def complement(sequence, reverse=False):
    """ Return complement of sequence """
    complements = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    sequence = reversed(sequence) if reverse else sequence
    return ''.join(complements.get(base, 'N') for base in sequence)
