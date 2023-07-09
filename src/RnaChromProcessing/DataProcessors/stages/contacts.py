from typing import Any, Dict, List

import numpy as np
import pandas as pd

from .basicstage import BasicStage
from ...utils import exit_with_error

DNA_COLUMNS = {
    'dna_chr': str, 'dna_bgn': np.uint32, 'dna_end': np.uint32,
    'id': str, 'dna_score': np.uint16, 'dna_strand': str, 'dna_cigar': str
}

RNA_COLUMNS = {
    'rna_chr': str, 'rna_bgn': np.uint32, 'rna_end': np.uint32,
    'id': str, 'rna_score': np.uint16, 'rna_strand': str, 'rna_cigar': str
}

MODES = ('fast', 'low-mem')

class Contacts(BasicStage):
    def __init__(self,
                 cfg: Dict[str, Any],
                 *args, **kwargs):
        super().__init__(*args, **kwargs)
        if (cpus := cfg.get('cpus', None)):
            self.cpus: int = cpus
        if ((mode := cfg.get('mode')) in MODES):
            self.mode: str = mode
        else:
            exit_with_error(f'Unknown operation mode for "contacts" stage: {mode=}')
        if self.mode == 'low-mem':
            self._chunksize: int = cfg.get('chunksize')

    def run(self,
            dna_ids: List[str],
            rna_ids: List[str]):
        """make contacts file from two bed files"""
        func = self._make_contacts if self.mode == 'fast' else self._make_contacts_chunks
        self.run_function(func, dna_ids, rna_ids)
    
    def _make_contacts(self,
                      dna_in_file: str,
                      rna_in_file: str,
                      _: str,  # ignore dna out file
                      rna_out_file: str) -> int:
        """read 2 bed files, produce 1 combined contacts file"""
         # prepare names for output and temporal files
        output = rna_out_file.rsplit('.', 1)[0] + '.tab'
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
    
    def _make_contacts_chunks(self,
                      dna_in_file: str,
                      rna_in_file: str,
                      _: str,  # ignore dna out file
                      rna_out_file: str) -> int:
        # https://stackoverflow.com/questions/58441517/merging-dataframe-chunks-in-pandas
        """read 2 bed files in a memory-efficient way, produce 1 combined contacts file"""
         # prepare names for output and temporal files
        output = rna_out_file.rsplit('.', 1)[0] + '.tab'
        dna = pd.read_csv(dna_in_file, sep='\t', header=None,
                          names=DNA_COLUMNS.keys(), dtype=DNA_COLUMNS)
        dna['id'] = dna['id'].apply(lambda x: x.split('.')[1] if '.' in x else x)
        result_chunks = []
        for chunk in pd.read_csv(rna_in_file, sep='\t', header=None,
                                 names=RNA_COLUMNS.keys(), dtype=RNA_COLUMNS,
                                 chunksize=self._chunksize):
            chunk['id'] = chunk['id'].apply(lambda x: x.split('.')[1] if '.' in x else x)
            result_chunks.append(pd.merge(dna, chunk, on='id', how='inner'))
        result = pd.concat(result_chunks)
        result = result.drop(['rna_score', 'dna_score'], axis=1)
        result.to_csv(output, sep='\t', index=False, header=True)
        return 0
