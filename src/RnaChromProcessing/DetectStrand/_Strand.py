import logging
import re
from os import chdir, listdir
from pathlib import Path
from tempfile import TemporaryDirectory, mkdtemp
from typing import Dict, List, Tuple

import pandas as pd
from pydantic import BaseModel, field_validator

from ..utils import check_file_exists, exit_with_error, find_in_list

BED_NAMES = ('chr', 'start', 'end', 'name', 'score', 'strand')
GTF_NAMES = ('chr', 'type', 'start', 'end', 'strand', 'attrs')

GENE_ID_PAT = re.compile(r'(?<=gene_id \")[^\"]+(?=\";)')
GENE_NAME_PAT = re.compile(r'(?<=gene_name \")[^\"]+(?=\";)')

logger = logging.getLogger('strand')
logging.getLogger("matplotlib").setLevel(logging.WARNING)


def _read_name_or_id(s: str) -> str:
    # will throw exception if both attrs ase absent (invalid annotation)
    capture = re.search(GENE_NAME_PAT, s) or re.search(GENE_ID_PAT, s)
    return capture.group(0)


class DetectStrand(BaseModel):
    prefix: str = 'strand'
    input_dir: Path
    output_dir: Path
    base_dir: Path = Path('.').resolve()

    gtf_annotation: Path
    genes_list: Path

    def __init__(self,
                 exp_groups: Dict[str, List[str]],
                 **kwargs):
        super().__init__(**kwargs)
        # directories
        self.__setup_dirs()
        chdir(self.base_dir)
        # tmp directory
        #self._work_dir = TemporaryDirectory(dir=self.base_dir)
        self._work_dir = mkdtemp(dir=self.base_dir)
        #self._work_pth = Path(self._work_dir.name)
        self._work_pth = Path(self._work_dir)
        # genes
        self.__load_genes()
        # input files
        self.__read_inputs(exp_groups)
        
    @field_validator('base_dir', 'input_dir', 'output_dir')
    @classmethod
    def resolve_path(cls, pth: str) -> Path:
        return Path(pth).resolve()
    
    @field_validator('gtf_annotation', 'genes_list')
    @classmethod
    def check_files(cls, pth: str) -> Path:
        check_file_exists(pth)
        return Path(pth).resolve()

    def __setup_dirs(self) -> None:
        if not self.bed_input_dir.exists():
            exit_with_error('Input directory does not exist!')
        if not self.base_dir.exists():
            self.base_dir.mkdir(parents=True)
        self.output_dir.mkdir(parents=True, exist_ok=True)
    
    def __load_genes(self) -> Path:
        # load lene names
        with open(self.genes_list, 'r') as f:
            # will search for gene_name "DDX11L2"; OR gene_id "ENSG00000290825.1"; for all listed genes
            genes_list: str = '|'.join([f'"{line.strip()}";' for line in f])
        # load annotation and extract selected genes
        gene_annot = pd.read_csv(
            self.gtf_annotation, sep='\t', header=None, skiprows=5,
            usecols=[0,2,3,4,6,8], names=GTF_NAMES
        )
        gene_annot = gene_annot[gene_annot['type'] == 'gene']
        gene_annot = gene_annot[gene_annot['attrs'].str.contains(genes_list)]
        if not (shape := gene_annot.shape[0]):
            exit_with_error('Could not find any of the listed genes in annotation file!')
        else:
            logger.info(f'{shape} genes selected from annotation file.')
        # save annotation as tmp bed file
        try:
            gene_annot['name'] = gene_annot['attrs'].apply(_read_name_or_id)
        except AttributeError:
            exit_with_error(
                'GTF annotation should have "gene_name" or "gene_id" attribute '
                ' for "gene" type records!'
            )
        gene_annot['score'] = 100
        self._bed_annot = self._work_pth / 'annotation.bed'
        gene_annot[BED_NAMES].to_csv(self._bed_annot, sep='\t', index=False, header=False)

    def __read_inputs(self, exp_groups: Dict[str, List[str]]) -> None:
        self._files_map: Dict[Tuple[str, str], str] = {}
        files_list = listdir(self.input_dir)
        for group, id_list in exp_groups.items():
            for file_id in id_list:
                filename = find_in_list(file_id, files_list)
                if not filename:
                    logger.warning(f'{file_id} not found in input directory, skipping..')
                    continue
                self._files_map[(group, file_id)] = filename
        if not self._files_map:
            exit_with_error('Could not find any if the listed files in input directory!')
        logger.info(f'{len(self._files_map)} files found in input directory')

    def run(self) -> None:
        # process files: tab or bed
        # for each file (in parallel)?:
        #   run bedtools coverage -s AND -S
        #   collect  and format results
        # organize results, make plots
        pass
