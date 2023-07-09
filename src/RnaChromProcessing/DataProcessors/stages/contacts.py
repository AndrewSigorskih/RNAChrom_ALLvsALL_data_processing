from typing import Any, Dict, List

import numpy as np
import pandas as pd

from .basicstage import BasicStage

DNA_COLUMNS = {
    'dna_chr': str, 'dna_bgn': np.uint32, 'dna_end': np.uint32,
    'id': str, 'dna_score': np.uint16, 'dna_strand': str, 'dna_cigar': str
}

RNA_COLUMNS = {
    'rna_chr': str, 'rna_bgn': np.uint32, 'rna_end': np.uint32,
    'id': str, 'rna_score': np.uint16, 'rna_strand': str, 'rna_cigar': str
}

DNA_COLUMN_NAMES = [
    'dna_chr', 'dna_bgn', 'dna_end', 'id',
    'dna_score', 'dna_strand', 'dna_cigar'
]

RNA_COLUMN_NAMES = [
    'rna_chr', 'rna_bgn', 'rna_end', 'id',
    'rna_score', 'rna_strand', 'rna_cigar'
]

class Contacts(BasicStage):
    def __init__(self,
                 cfg: Dict[str, Any],
                 *args, **kwargs):
        super().__init__(*args, **kwargs)
        if (cpus := cfg.get('cpus', None)):
            self.cpus: int = cpus

    def run(self,
            dna_ids: List[str],
            rna_ids: List[str]):
        """make contacts file from two bed files"""
        func = self._make_contacts
        self.run_function(func, dna_ids, rna_ids)
    
    def _make_contacts(self,
                      dna_in_file: str,
                      rna_in_file: str,
                      _: str,  # ignore dna out file
                      rna_out_file: str):
        """read 2 bed files, produce 1 combined contacts file"""
         # prepare names for putput and temporal files
        output = rna_out_file.rsplit('.', 1)[0]+ '.tab'
        dna = pd.read_csv(dna_in_file, sep='\t', header=None,
                          names=DNA_COLUMNS.keys(), dtype=DNA_COLUMNS)
        rna = pd.read_csv(rna_in_file, sep='\t', header=None,
                          names=RNA_COLUMNS.keys(), dtype=RNA_COLUMNS)
        dna['id'] = dna['id'].apply(lambda x: x.split('.')[1] if '.' in x else x)
        rna['id'] = rna['id'].apply(lambda x: x.split('.')[1] if '.' in x else x)
        rna = pd.merge(rna, dna, on='id', how='inner')
        rna = rna.drop(['rna_score', 'dna_score'], axis=1)
        rna.to_csv(output, sep='\t', index=False, header=True)
        return 0
