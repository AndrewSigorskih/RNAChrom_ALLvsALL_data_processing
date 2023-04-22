import logging
import os

import pandas as pd

CONTACTS_COLS = ('rna_chr', 'rna_bgn', 'rna_end', 'rna_strand')
GTF_NAMES = ('chr', 'type', 'bgn', 'end', 'strand', 'attrs')
logger = logging.getLogger()

class StrandCalc:
    def __init__(self,
                 input_dir: str,
                 gtf_annotation: str,
                 genes: str,
                 output_dir: str,
                 prefix: str) -> None:
        self.input_dir = input_dir  # not needed?
        self.output_dir = output_dir
        self.prefix = prefix
        self.files = [
            os.path.join(input_dir, file) for file in os.listdir(input_dir)
            if file.endswith('.tab')  # subject to change based on mapping provided
        ]
        self._load_genes(genes, gtf_annotation)
    
    def _load_genes(self,
                    genes: str,
                    gtf_annotation: str) -> None:
        with open(genes, 'r') as f:
            # will search for gene_name "DDX11L2"; OR gene_id "ENSG00000290825.1"; for all listed genes
            genes_list: str = '|'.join([f'"{line.strip()}";' for line in f])
        gene_annot = pd.read_csv(
            gtf_annotation, sep='\t', header=None, skiprows=5,
            usecols=[0,2,3,4,6,8], names=GTF_NAMES
        )
        gene_annot = gene_annot[gene_annot['type'] == 'gene']
        gene_annot = gene_annot[gene_annot['attrs'].str.contains(genes_list)]
        logger.debug(f'{gene_annot.shape[0]} genes selected from annotation file.')
        gene_annot['idx'] = gene_annot['attrs'].apply(lambda x: x.split('gene_name')[1].split(';')[0].split('"')[1])
        gene_annot = gene_annot.set_index('idx')
        self.gene_annot: pd.DataFrame = gene_annot

    def calculate(self) -> None:
        result = pd.DataFrame(data=None, columns=self.gene_annot.index, index=self.files)
        for inp_file in self.files:
            logger.debug(f'started processing {inp_file}')
            data = pd.read_csv(inp_file, sep='\t', usecols=CONTACTS_COLS)
            for gene in result.columns:
                mask = ((data['rna_chr'] == self.gene_annot.at[gene, 'chr']) &
                        (data['rna_bgn'] <= self.gene_annot.at[gene, 'end']) &
                        (data['rna_end'] >= self.gene_annot.at[gene, 'bgn'])
                )
                same = (mask & (data['rna_strand'] == self.gene_annot.at[gene, 'strand'])).sum()
                anti = (mask & (data['rna_strand'] != self.gene_annot.at[gene, 'strand'])).sum()
                result.at[inp_file, gene] = (same, anti)
        self.result: pd.DataFrame = result


    def run(self) -> None:
        # https://stackoverflow.com/questions/17322109/get-dataframe-row-count-based-on-conditions
        # calc results and store in table
        # get json for further xrna properties
        # make plots
        # save outputs
        self.calculate()
        self.result.to_csv(f'{self.output_dir}/{self.prefix}.tsv', sep='\t')
    '''
    for name in names:
        dat = pd.read_csv(f'~/data/imargi/bed/{name}.rna.bed', sep='\t', header=None)
        print(name, dat.shape[0])
        dat[0] = dat[0].apply(lambda x: hg38_dict[x])
        for item in ribgenes.index:
            same = dat[(dat[0]==ribgenes.at[item,0]) & (dat[5]==ribgenes.at[item,6]) & (dat[1] < ribgenes.at[item,4]) & (dat[2] > ribgenes.at[item,3])].shape[0]
            anti = dat[(dat[0]==ribgenes.at[item,0]) & (dat[5]!=ribgenes.at[item,6]) & (dat[1] < ribgenes.at[item,4]) & (dat[2] > ribgenes.at[item,3])].shape[0]
            res.at['imargi_'+name,item] = [same, anti]
    '''
    # save several otputs
    # save images (png and pdf)