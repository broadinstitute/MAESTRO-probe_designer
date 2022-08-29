import click
from copy import copy
import logging
import pandas as pd
import numpy as np

from probe_designer import __version__
from .sequence import Probe
from .utils import FastaReader
from .thermodynamics import ThermoParams

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)


logging.basicConfig(
    format='%(levelname)s\t%(asctime)s\t%(name)s\t%(message)s',
    level=logging.INFO)
logger = logging.getLogger('ProbeDesigner')


def maf_to_df(maf_file):
    """ Load MAF, return dataframe """
    df = pd.read_table(maf_file, dtype={'Chromosome': str}, comment="#")
    df.columns = df.columns.str.lower()
    logger.info(f'Found {len(df)} targets in MAF file')
    return df


def sort_maf(df):
    """ Sort chromosome, start positions correctly (1-22, X, Y, etc.) """
    df['chromosome'] = df['chromosome'].apply(
        lambda x: int(x) if x.isdigit() else x)
    df = df.sort_values(['chromosome', 'start_position'])
    df['chromosome'] = df['chromosome'].astype(str)
    df.reset_index(inplace=True, drop=True)
    return df


def remove_duplicates(df):
    """ Remove duplicates from MAF dataframe """
    with_duplicates = len(df)
    check_dup_cols = ['chromosome', 'start_position', 'tumor_seq_allele2']
    df.drop_duplicates(subset=check_dup_cols, inplace=True)
    duplicate_counts = with_duplicates - len(df)
    if duplicate_counts > 0:
        logger.info(f'Removed {duplicate_counts} duplicates from MAF')


def process_maf(maf_file, snps_only=True):
    """ Load MAF, remove duplicates and remove non-SNPs, and return
    iterator of MAF entries as a namedtuple """
    maf_df = maf_to_df(maf_file)
    remove_duplicates(maf_df)
    if snps_only and 'variant_type' in maf_df.columns:
        total = len(maf_df)
        maf_df = maf_df[maf_df.variant_type == 'SNP']
        not_snp = total - len(maf_df)
        if not_snp > 0:
            logger.info(f'Removed {not_snp} targets that were not SNPs')

    maf_df = sort_maf(maf_df)
    return maf_df


def create_probe(variant, probe_length, fasta, thermo, dg_range,
                 distance, static_length, design_ref, same_strand,
                 disable_blast):
    """ Create probe from variant object.
    There are 2 strategies for creating a probe:
    Static length - all probes have same length with variant in middle of probe
    Variable length - optimal probe is selected based on given delta G range
    """
    probe = Probe(
        chrom=variant.chromosome,
        alt_position=variant.start_position,
        ref_base=variant.reference_allele,
        alt_base=variant.tumor_seq_allele2,
        length=probe_length,
        reference=fasta,
        thermo=thermo,
        design_ref=design_ref,
        same_strand=same_strand)
    if static_length:
        if not disable_blast:
            probe.get_filtered_blast()
        return probe
    # Slide a window across target base, and grab longest probe for each
    # window. By designing the longest probes, we should be using the
    # most unique sequence (thus least BLAST hits)
    # We try only 3 conditions here to limit the amount of BLAST
    # searches we need to perform for each target
    options = []
    for i in [0, -distance, distance]:
        adj_probe = copy(probe)
        adj_probe.set_offset(i)
        adj_probe = get_longest_length(adj_probe, dg_range)
        options.append(adj_probe)
    return get_min_blast(options, disable_blast)


def create_probes(maf_df, probe_length, fasta, thermo, dg_range,
                  distance, static_length, design_ref, same_strand,
                  disable_blast):
    """ Create probes from MAF """
    logger.info('Creating probes and running BLAST...')
    for _, variant in maf_df.iterrows():
        yield create_probe(variant, probe_length, fasta, thermo, dg_range,
                           distance, static_length, design_ref, same_strand,
                           disable_blast)


def get_longest_length(probe, dg_range):
    template = copy(probe)
    dg_min, dg_max = dg_range
    # non-ACGT bases cause dg = nan
    if np.isnan(template.dg):
        return template
    if template.dg < dg_min:
        while template.dg < dg_min:
            previous = copy(template)
            # Edge case with large offsets and shortening length
            if abs(template.offset) >= 0.5 * (template.length - 1):
                longest_probe = previous
                break
            template.set_length(template.length - 1)
            longest_probe = copy(template)
    else:
        while template.dg >= dg_min:
            previous = copy(template)
            longest_probe = previous
            template.set_length(template.length + 1)
    if (longest_probe.dg < dg_min) or (longest_probe.dg > dg_max):
        # When designing, if a probe cannot be designed within
        # dg_range, probe closest to dg_max will be chosen
        options = [longest_probe, previous]
        longest_probe = min(options, key=lambda x: abs(x.dg - dg_min))
    return longest_probe


def get_min_blast(options, disable_blast):
    """ List of probe objects """
    if not disable_blast:
        for probe in options:
            probe.get_filtered_blast()
    best_probe = min(options, key=lambda x: (x.blast_count, abs(x.offset)))
    return best_probe


def probe_to_string(probe, attributes):
    info = map(str, [getattr(probe, attr) for attr in attributes])
    return '\t'.join(info)


class FloatRange(click.ParamType):
    name = 'floatrange'

    def convert(self, value, param, ctx):
        try:
            return list(map(float, value.split(',')))
        except ValueError:
            self.fail(
                f'{value} is not a comma separated float range',
                param,
                ctx,
            )


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option(
    '-m',
    '--maf',
    type=click.Path(exists=True),
    help='MAF file for design',
    required=True)
@click.option(
    '-o',
    '--output',
    type=click.Path(writable=True),
    help='Output file',
    required=True)
@click.option(
    '-r',
    '--reference',
    type=click.Path(exists=True),
    help='Reference file [Broad HG19]',
    default='/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta')
@click.option(
    '-p',
    '--probe_length',
    type=click.INT,
    metavar='INT',
    help='Probe length in bp',
    default=30,
    show_default=True)
@click.option(
    '-t',
    '--temp',
    help='Temperature of annealing',
    type=click.FLOAT,
    default=50.0,
    show_default=True)
@click.option(
    '-d',
    '--dg_range',
    help='Min, max delta G of probe for design',
    type=FloatRange(),
    default='-18.0,-14.0',
    show_default=True)
@click.option(
    '--static_length',
    help='Design static length probe (ignores delta G)',
    is_flag=True
)
@click.option(
    '--same_strand',
    help='Use same strand for all probes',
    is_flag=True
)
@click.option(
    '--design_ref',
    help='Design probe with reference base',
    is_flag=True
)
@click.option(
    '--disable_blast',
    help='Run probe design without running BLAST',
    is_flag=True
)
@click.version_option(None, "-v", "--version", message="%(version)s")
def probe_designer(maf, output, reference, probe_length, temp, dg_range,
                   static_length, same_strand, design_ref, disable_blast):
    """ Probe designer for MAESTRO """
    arguments = '\n'.join([f'{k}={v}' for k, v in locals().items()])
    version = __version__
    logger.info(f'Running using {version} with arguments:\n\n{arguments}\n')
    thermo = ThermoParams(temp=temp, Na=50.0, Mg=0.0, dnac1=250.0, dnac2=250.0)
    fasta = FastaReader(reference)
    maf_df = process_maf(maf)
    distance = 5
    probes = create_probes(maf_df, probe_length, fasta,
                           thermo, dg_range, distance, static_length,
                           design_ref, same_strand, disable_blast)
    # All header names must be probe attributes
    header = [
        'chrom', 'start', 'end', 'alt_position', 'ref_base', 'alt_base',
        'sequence', 'strand', 'length', 'gc', 'tm', 'dg', 'blast_count'
    ]
    with open(output, 'w') as outfile:
        outfile.write('\t'.join(header) + '\n')
        count = 0
        for probe in probes:
            probe_string = probe_to_string(probe, header)
            count += 1
            outfile.write(probe_string + '\n')
            logger.info(f'Completed design for {count} probes')


if __name__ == '__main__':
    probe_designer()
