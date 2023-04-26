import json
import logging
import os
from typing import Dict, List

import pandas as pd

from ..utils import check_file_exists, exit_with_error, find_in_list, make_directory
from ..plots import rna_strand_barplot, set_style_white

CONFIG_FIELDS = ('input_dir', 'output_dir', 'gtf_annotation', 'genes_list', 'exp_groups')
CONTACTS_COLS = ('rna_chr', 'rna_bgn', 'rna_end', 'rna_strand')
GTF_NAMES = ('chr', 'type', 'bgn', 'end', 'strand', 'attrs')

logger = logging.getLogger('strand')
logging.getLogger("matplotlib").setLevel(logging.WARNING)

class StrandCalc:
    def __init__(self,
                 config_file: str) -> None:
        
        with open(config_file, 'r') as f:
            config: dict = json.load(f)

        # check config contents
        missing = []
        for field in CONFIG_FIELDS:
            if field not in config.keys():
                missing.append(field)
        if missing:
            exit_with_error(f'Required fields not found in config: {", ".join(missing)}')
        
        # check files exist and load annotation data
        check_file_exists(config['genes_list'])
        check_file_exists(config['gtf_annotation'])
        self._load_genes(config['genes_list'], config['gtf_annotation'])

        # set variables
        self.prefix: str = config.get('prefix', 'strand')
        self.input_dir: str = config['input_dir']
        self.output_dir: str = config['output_dir']
        make_directory(self.output_dir)

        # check for files
        self.files_map: Dict[str, str] = {}
        exp_groups: Dict[str, List[str]] = config['exp_groups']
        files_list = os.listdir(self.input_dir)
        for group, id_list in exp_groups.items():
            for file_id in id_list:
                filename = find_in_list(file_id, files_list)
                if not filename:
                    logger.warning(f'{file_id} not found in input directory, skipping..')
                    continue
                self.files_map[f'{group}_{file_id}'] = filename
        if not self.files_map:
            exit_with_error('Could not find any if the listed files in input directory!')
        logger.info(f'{len(self.files_map)} files found in input directory')
    
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
        if not (shape := gene_annot.shape[0]):
            exit_with_error('Could not find any of the listed genes in annotation file!')
        else:
            logger.info(f'{shape} genes selected from annotation file.')
        gene_annot['idx'] = gene_annot['attrs'].apply(lambda x: x.split('gene_name')[1].split(';')[0].split('"')[1])
        gene_annot = gene_annot.set_index('idx')
        self.gene_annot: pd.DataFrame = gene_annot

    def calculate(self) -> None:
        result = pd.DataFrame(data=None,
                              columns=self.gene_annot.index,
                              index=self.files_map.keys())
        for name, file in self.files_map.items():
            logger.debug(f'Started processing {file}')
            data = pd.read_csv(f'{self.input_dir}/{file}', sep='\t', usecols=CONTACTS_COLS)
            for gene in result.columns:
                mask = ((data['rna_chr'] == self.gene_annot.at[gene, 'chr']) &
                        (data['rna_bgn'] <= self.gene_annot.at[gene, 'end']) &
                        (data['rna_end'] >= self.gene_annot.at[gene, 'bgn'])
                )
                same: int = (mask & (data['rna_strand'] == self.gene_annot.at[gene, 'strand'])).sum()
                anti: int = (mask & (data['rna_strand'] != self.gene_annot.at[gene, 'strand'])).sum()
                result.at[name, gene] = (same, anti)
        self.raw_result: pd.DataFrame = result

    def counts_to_values(self) -> None:
        same_wins = lambda x: x[0] > x[1]
        anti_wins = lambda x: x[1] > x[0]
        self.result = pd.DataFrame(data=None, index=self.raw_result.index)
        self.result['same'] = self.raw_result.apply(lambda row: row.apply(same_wins).sum(),
                                                    axis=1)
        self.result['anti'] = self.raw_result.apply(lambda row: row.apply(anti_wins).sum(),
                                                    axis=1)
        
        mask = (self.result['same'] - self.result['anti'])/(self.result['same'] + self.result['anti'])
        self.result['strand'] = 'UNKNOWN'
        self.result.loc[mask > 0.75, 'strand'] = 'SAME'
        self.result.loc[mask < -0.75, 'strand'] = 'ANTI'


    def run(self) -> None:
        # calc results and store in table
        self.calculate()
        self.counts_to_values()
        # save outputs
        self.result.to_csv(f'{self.output_dir}/{self.prefix}_wins.tsv', sep='\t')
        self.raw_result.to_csv(f'{self.output_dir}/{self.prefix}_raw_counts.tsv', sep='\t')
        # make plots
        set_style_white()
        rna_strand_barplot(self.result, self.gene_annot.shape[0],
                           self.output_dir, self.prefix)
        logger.info('Done.')
