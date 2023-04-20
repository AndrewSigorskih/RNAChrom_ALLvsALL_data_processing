import logging
from typing import List

import pandas as pd

GTF_NAMES = ('chr', 'type', 'start', 'end', 'strand', 'attrs')
logger = logging.getLogger()

class StrandCalc:
    def __init__(self,
                 input_dir: str,
                 gtf_annotation: str,
                 genes: str,
                 output_dir: str,
                 prefix: str) -> None:
        self.input_dir = input_dir
        self.output_dir = output_dir
        self.prefix = prefix
        self._load_genes(genes, gtf_annotation)
    
    def _load_genes(self,
                    genes: str,
                    gtf_annotation: str) -> None:
        with open(genes, 'r') as f:
            genes_list: List[str] = [line.strip() for line in f]
            self.gene_annot: pd.DataFrame = pd.read_csv(
                gtf_annotation, sep='\t', header=None, skiprows=5,
                usecols=[0,2,3,4,6,8], names=GTF_NAMES
            )
            self.gene_annot = self.gene_annot[self.gene_annot['type'] == 'gene']

    def run(self) -> None:
        pass