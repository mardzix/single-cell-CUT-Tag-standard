#!/usr/bin/env python3

import argparse
import snakemake

parser = argparse.ArgumentParser(description='Preprocessing pipeline for single-cell CUT&Tag')

parser.add_argument('-i', '--cellranger',
                    dest='cellranger_path',
                    required=True,
                    type=str,
                    help='path to cellranger output folder')

parser.add_argument('-o', '--out',
                    dest='out_prefix',
                    required=True,
                    type=str,
                    help='path to cellranger output folder')



parser.add_argument('-s', '--sample_name',
                    required=True,
                    type=str,
                    help='Sample name name')

parser.add_argument('-m,', '--min_reads_log10',
                    type=float,
                    default=4.0,
                    help='minimum number of reads per cells to pass the filter [default 4.0]')

parser.add_argument('-n', '--max_reads_log10',
                    type=float,
                    default=5.5,
                    help='maximum number of reads per cells to pass the filter [default 5.5]')

parser.add_argument('-f', '--min_frip_value',
                    type=float,
                    default=0.25,
                    help='minimum fraction of reads in peaks to pass the filter [default 0.25]')

parser.add_argument('-g', '--max_frip_value',
                    type=float,
                    default=0.95,
                    help='maximum fraction of reads in peaks to pass the filter [default 0.95]')

parser.add_argument('-a', '--antibody',
                    type=str,
                    default='Unknown_antibody',
                    help='Antibody name')



parser.add_argument('-b', '--bin_window',
                    type=int,
                    default=5000,
                    help='Bin width to use for matrix construction [default 5000]')

parser.add_argument('-x', '--genome_version',
                    type=str,
                    default='mm10',
                    help='Version of the genome [currently only mm10 is supported; default mm10]')

args = parser.parse_args()
print(args)

config = vars(args)

snakemake.snakemake(snakefile = 'workflow/snakefile_preprocess',
                    cores=4,
                    config=config)
