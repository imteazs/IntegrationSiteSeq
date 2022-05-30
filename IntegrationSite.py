import pandas as pd
from Bio import SeqIO
import argparse
from pathlib import Path


class IntegrationSite:
    def __init__(self, gtfdf):
        self.gtfdf = gtfdf

    'pull the sequence with the given start, stop, chr, and genome dictionary'
    def pullSeq(self, start, stop, chr_name, genedict):
        record = genedict[chr_name]
        seq = str(record.seq)[start:stop]
        return seq


if __name__ == '__main__' :
    parser = argparse.ArgumentParser()
    parser.add_argument('--gtf_path', help='Path to gtf file, this location is where the output csv is saved',
                        required=True)
    parser.add_argument('--genome_path', help='Path to genome', required=True)

    args = parser.parse_args()

    'output path'
    output = Path(args.gtf_path).parent.joinpath('result_seq.csv')

    'load a temp dataframe and rename columns'
    integratedf = pd.read_csv(args.gtf_path, sep='\t', comment='#', header=None)
    gtf_columns = ['chr_name', 'source', 'feature_type', 'start', 'stop', 'score', 'strand', 'genom_phase', 'add_info']
    integratedf.columns = gtf_columns

    'load genome as a dictionary to efficiently parse through'
    genome_dict = SeqIO.to_dict(SeqIO.parse(args.genome_path, 'fasta'))

    'Massaging the dataframe'
    integratedf = integratedf.loc[integratedf['feature_type'] == 'stop_codon']
    integratedf = integratedf[['chr_name', 'source', 'feature_type', 'start']].reset_index(drop=True)
    chr_start = integratedf['start'].iloc[1:].reset_index(drop=True)
    integratedf = integratedf[:-1]
    integratedf = integratedf.join(chr_start, rsuffix='_next_stop_codon')

    'Shift start poisition back by one, due to python indexing starting at 0'
    integratedf['pystart_stop_codon'] = integratedf['start'] - 1
    integratedf['pystart_next_stop_codon'] = integratedf['start_next_stop_codon'] - 1

    selectdf = IntegrationSite(integratedf)
    'Use the pullSeq function to get the sequence'
    selectdf.gtfdf['sequence'] = selectdf.gtfdf.apply(lambda x: selectdf.pullSeq(x['pystart_stop_codon'],
                                                                                 x['pystart_next_stop_codon'],
                                                                                 x['chr_name'], genome_dict), axis=1)

    'remove null values'
    selectdf.gtfdf = selectdf.gtfdf.loc[selectdf.gtfdf['sequence'] != '']
    'Getting sequence length'
    selectdf.gtfdf['seq_len'] = selectdf.gtfdf['sequence'].apply(len)
    'Check whether sequences has lowercase nucleotides'
    selectdf.gtfdf['has_lowercase'] = selectdf.gtfdf['sequence'].str.contains('[atcg]')
    'Check if string starts with stop codon TAG, TAA, TGA'
    selectdf.gtfdf['seq_starts_tag_taa_tga'] = selectdf.gtfdf['sequence'].str.contains('^(tag)|(taa)|(tga)')

    'save file'
    selectdf.gtfdf.to_csv(output, index=False)
