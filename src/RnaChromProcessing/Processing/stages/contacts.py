from dataclasses import dataclass
from pathlib import Path
from typing import List, Literal

import numpy as np
import pandas as pd
from pydantic import PositiveInt

from .basicstage import BasicStage, SamplePair
from ...utils import exit_with_error

DNA_COLUMNS = {
    'dna_chr': str, 'dna_bgn': np.uint32, 'dna_end': np.uint32,
    'id': str, 'dna_score': np.uint16, 'dna_strand': str, 'dna_cigar': str
}

RNA_COLUMNS = {
    'rna_chr': str, 'rna_bgn': np.uint32, 'rna_end': np.uint32,
    'id': str, 'rna_score': np.uint16, 'rna_strand': str, 'rna_cigar': str
}

TAB_HEADER = (
    'rna_chr', 'rna_bgn', 'rna_end', 'id', 'rna_strand', 'rna_cigar',
    'dna_chr', 'dna_bgn', 'dna_end', 'dna_strand', 'dna_cigar'
)


def _process_id(read_id: str) -> str:
    """remove 'SRRXXXXX' part from read ID if exists"""
    if '.' in read_id:
        return read_id.split('.', maxsplit=1)[1]
    return read_id


@dataclass
class BedRow:
    chr: str
    bgn: str
    end: str
    id: str
    score: str
    strand: str
    cigar: str

    @classmethod
    def from_string(cls, row: str):
        return cls(*row.split())._prepare_id()

    def _prepare_id(self):
        self.id = _process_id(self.id)
        return self

    def to_rna(self) -> str:
        return '\t'.join(
            self.chr, self.bgn, self.end, self.id, self.strand, self.cigar
        )
    
    def to_dna(self) -> str:
        return '\t'.join(
            self.chr, self.bgn, self.end, self.strand, self.cigar
        )


class Contacts(BasicStage):
    mode: Literal['fast', 'chunks', 'low-mem'] = 'low-mem'
    chunksize: PositiveInt = 10_000_000

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
    
    def run(self, samples: List[SamplePair]) -> None:
        """make tab-separated contacts file from two bed files"""
        # choose function to run
        if self.mode == 'fast':
            func = self._make_contacts_fast
        elif self.mode == 'chunks':
            func = self._make_contacts_chunks
        elif self.mode == 'low-mem':
            func = self._make_contacts_iter
        # prepare filepaths
        output_samples = self._make_output_samples(samples)
        # run function
        self.run_function(
            func,
            [sample.dna_file for sample in samples],
            [sample.rna_file for sample in samples],
            [sample.dna_file for sample in output_samples],
            [sample.rna_file for sample in output_samples]
        )

    def _make_contacts_fast(self,
                            dna_in_file: Path,
                            rna_in_file: Path,
                            _: Path,  # ignore dna out file
                            rna_out_file: Path) -> int:
        """read 2 bed files in memory, produce 1 combined contacts file"""
        dna = pd.read_csv(dna_in_file, sep='\t', header=None,
                          names=DNA_COLUMNS.keys(), dtype=DNA_COLUMNS)
        rna = pd.read_csv(rna_in_file, sep='\t', header=None,
                          names=RNA_COLUMNS.keys(), dtype=RNA_COLUMNS)
        dna['id'] = dna['id'].apply(_process_id)
        rna['id'] = rna['id'].apply(_process_id)
        rna = pd.merge(rna, dna, on='id', how='inner')
        rna = rna.drop(['rna_score', 'dna_score'], axis=1)
        rna.to_csv(rna_out_file, sep='\t', index=False, header=True)
        return 0
    
    def _make_contacts_chunks(self,
                              dna_in_file: Path,
                              rna_in_file: Path,
                              _: Path,  # ignore dna out file
                              rna_out_file: Path) -> int:
        # https://stackoverflow.com/questions/58441517/merging-dataframe-chunks-in-pandas
        """read rna bed file in a memory-efficient way, produce 1 combined contacts file"""
         # prepare names for output and temporal files
        dna = pd.read_csv(dna_in_file, sep='\t', header=None,
                          names=DNA_COLUMNS.keys(), dtype=DNA_COLUMNS)
        dna['id'] = dna['id'].apply(_process_id)
        result_chunks = []
        for chunk in pd.read_csv(rna_in_file, sep='\t', header=None,
                                 names=RNA_COLUMNS.keys(), dtype=RNA_COLUMNS,
                                 chunksize=self.chunksize):
            chunk['id'] = chunk['id'].apply(_process_id)
            result_chunks.append(pd.merge(dna, chunk, on='id', how='inner'))
        result = pd.concat(result_chunks)
        result = result.drop(['rna_score', 'dna_score'], axis=1)
        result.to_csv(rna_out_file, sep='\t', index=False, header=True)
        return 0
    
    def _make_contacts_iter(self,
                            dna_in_file: Path,
                            rna_in_file: Path,
                            _: Path,  # ignore dna out file
                            rna_out_file: Path) -> int:
        """read 2 bed files line-by-line, produce 1 combined contacts file"""
        with open(dna_in_file, 'r') as dna_in, \
             open(rna_in_file, 'r') as rna_in, \
             open(rna_out_file, 'w') as rna_out:
            print(*TAB_HEADER, file=rna_out, sep='\t')
            rna_row, dna_row = next(rna_in), next(dna_in)
            while(rna_row and dna_row):
                rna_obj, dna_obj = BedRow.from_string(rna_row), BedRow.from_string(dna_row)
                while rna_obj.id < dna_obj.id:
                    rna_obj = BedRow.from_string(next(rna_in).strip())
                while rna_obj.id > dna_obj.id:
                    dna_obj = BedRow.from_string(next(dna_in).strip())
                print(rna_obj.to_rna(), dna_obj.to_dna(), file=rna_out, sep='\t')
                rna_row, dna_row = next(rna_in), next(dna_in)
        return 0
