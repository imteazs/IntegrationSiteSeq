import pandas as pd
from Bio import SeqIO
import argparse
from pathlib import Path


class IntegrationSite:
    def __init__(self, gtfdf):
        self.gtfdf = gtfdf


if __name__ == '__main__' :
    parser = argparse.ArgumentParser()
    parser.add_argument('--gtf_path', help='Path to gtf file, this location is where the output csv is saved', required=True)
    parser.add_argument('--genome_path', help='Path to genome', required=True)

    args = parser.parse_args()

    output = Path(args.gtf_path).parent.joinpath('result_seq.csv')

    integratedf = pd.read_csv(args.gtf_path, sep='\t', comment='#', header=None)
    gtf_columns = ['chr_name', 'source', 'feature_type', 'start', 'stop', 'score', 'strand', 'genom_phase', 'add_info']
    integratedf.columns = gtf_columns
