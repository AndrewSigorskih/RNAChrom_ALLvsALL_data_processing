from typing import List

import pandas as pd

from .basicstage import BasicStage

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
                 *args, **kwargs):
        super().__init__(*args, **kwargs)

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
                          names=DNA_COLUMN_NAMES)
        rna = pd.read_csv(rna_in_file, sep='\t', header=None,
                          names=RNA_COLUMN_NAMES)
        dna['id'] = dna['id'].apply(lambda x: x.split('.')[1] if '.' in x else x)
        rna['id'] = rna['id'].apply(lambda x: x.split('.')[1] if '.' in x else x)
        rna = pd.merge(rna, dna, on='id', how='inner')
        rna = rna.drop(['rna_score', 'dna_score'], axis=1)
        rna.to_csv(output, sep='\t', index=False, header=True)
        return 0
