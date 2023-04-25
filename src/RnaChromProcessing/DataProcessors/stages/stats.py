import concurrent.futures
import logging
import os
from itertools import repeat
from typing import Callable, List, Optional

import pandas as pd

from .basicstage import suff_to_filter
from ...utils import find_in_list, run_command, run_get_stdout

logger = logging.getLogger()
how_variants = ('skip', 'default', 'full')

class StatsCalc:
    def __init__(self, 
                 output_dir: str,
                 cpus: int,
                 how: str,
                 prefix: str,
                 dna_ids: List[str],
                 rna_ids: List[str]):
        self.output_dir: str = output_dir
        self.cpus: int = cpus
        self.dna_ids = dna_ids
        self.rna_ids = rna_ids
        self.prefix = prefix
        self.result = {}
        os.makedirs('stats', exist_ok=True)
        if how not in how_variants:
            logger.warning(f'Unknown option for read statistics calcultaion: {how}. Switching to default.')
            self.how = 'default'
        else:
            self.how = how

    def run_function(self,
                     func: Callable[[str, str, Optional[str]], int],
                     folder: str,
                     rna_input_files: List[str],
                     dna_input_files: List[str],
                     mode: Optional[str] = None):
        logger.debug(f'Calculating surviving reads statistic in {folder}')
        key = mode if mode else folder
        with concurrent.futures.ThreadPoolExecutor(max_workers=self.cpus) as executor:
            future_to_id = {
                executor.submit(func, rna_file, dna_file, mode) : rna_id
                for rna_file, dna_file, rna_id in zip(rna_input_files, dna_input_files, self.rna_ids)
            }
            for future in concurrent.futures.as_completed(future_to_id):
                rna_id = future_to_id[future]
                count = future.result()
                self.result[key][rna_id] = count

    def count_in_fastq_pair(self,
                            rna_file: str,
                            dna_file: str,
                            *_) -> int:
        """Counts matching IDs in pair of FASTQ files.
        Supported IDs formats:
            * @IDXXX <-> @IDXXX\n
            * @FILE1.IDXXX <-> @FILE2.IDXXX\n
            * @FILE.IDXXX <-> @FILE.IDXXX"""
        dna_tmp_file = os.path.join('stats', rna_file.split('/')[-1] + '.tmp')
        rna_tmp_file = os.path.join('stats', dna_file.split('/')[-1] + '.tmp')

        for infile, outfile in ((rna_file, rna_tmp_file), (dna_file, dna_tmp_file)):
            cat = 'zcat' if infile.endswith('.gz') else 'cat'
            cmd = (
                f'{cat} {infile} | sed -n  "1~4p" | sed "s/@//g"  | '
                f'cut -d " " -f 1 | cut -d "." -f 2 > {outfile}'
            )
            run_command(cmd, shell=True)

        count_cmd = f"awk 'a[$0]++' {dna_tmp_file} {rna_tmp_file} | wc -l"
        result = run_get_stdout(count_cmd, shell=True)
        return int(result)
    
    def count_in_bam_pair(self,
                          rna_file: str,
                          dna_file: str,
                          mode='mapped') -> int:
        """Counts matching IDs in pair od BAM files.
        Supported IDs formats:
            * @IDXXX <-> @IDXXX\n
            * @FILE1.IDXXX <-> @FILE2.IDXXX\n
            * @FILE.IDXXX <-> @FILE.IDXXX"""
        dna_tmp_file = os.path.join('stats', rna_file.split('/')[-1] + '.tmp')
        rna_tmp_file = os.path.join('stats', dna_file.split('/')[-1] + '.tmp')
        for infile, outfile in ((rna_file, rna_tmp_file), (dna_file, dna_tmp_file)):
            if mode == 'mapped':
                cmd = (
                    f"samtools view -F 4 {infile} | "
                    "awk -F \"\\t\" '{print \"@\"$1}' | "
                    #"cut -d '.' -f2 | awk '{ a[$1]++ } END { for (b in a) { print b } }' " # change to "sort -u"
                    "cut -d '.' -f2 | sort -u "
                    f"> {outfile}"
                )
            elif mode == 'mapped2mism':
                cmd = (
                    f"samtools view -F 4 {infile} | grep -E 'XM:i:[0-2]\s.*' | "
                    "awk -F \"\\t\" '{print \"@\"$1}' | cut -d '.' -f2 | "
                    #"awk '{ a[$1]++ } END { for (b in a) { print b } }' "
                    "sort -u "
                    f'> {outfile}'
                )
            else:
                raise ValueError(f'Unknown mode for count_in_bam_pair: {mode}')
            run_command(cmd, shell=True)
        count_cmd = f"awk 'a[$0]++' {dna_tmp_file} {rna_tmp_file} | wc -l"
        result = run_get_stdout(count_cmd, shell=True)
        return int(result)

    def count_in_fastqs(self, folder: str) -> None:
        self.result[folder] = {}
        # get and prepare filenames
        filenames: List[str] = [x for x in os.listdir(folder)
                                if not any([x.endswith(y) for y in suff_to_filter])]
        dna_files = [find_in_list(id, filenames) for id in self.dna_ids]
        rna_files = [find_in_list(id, filenames) for id in self.rna_ids]
        dna_input_files = [os.path.join(folder, filename)
                           for filename in dna_files]
        rna_input_files = [os.path.join(folder, filename)
                           for filename in rna_files]
        # concurrently run and fill result dict
        self.run_function(self.count_in_fastq_pair, folder, 
                          rna_input_files, dna_input_files)
    
    def count_in_bams(self, folder: str, mode: str) -> None:
        self.result[mode] = {}
        filenames: List[str] = [x for x in os.listdir(folder)
                                if not any([x.endswith(y) for y in suff_to_filter])]
        dna_files = [find_in_list(id, filenames) for id in self.dna_ids]
        rna_files = [find_in_list(id, filenames) for id in self.rna_ids]
        dna_input_files = [os.path.join(folder, filename)
                           for filename in dna_files]
        rna_input_files = [os.path.join(folder, filename)
                           for filename in rna_files]
        self.run_function(self.count_in_bam_pair, folder, 
                          rna_input_files, dna_input_files, mode)
        
    def count_contacts(self, folder: str) -> None:
        self.result[folder] = {}
        filenames: List[str] = [x for x in os.listdir(folder)]
        rna_files = (find_in_list(id, filenames) for id in self.rna_ids)
        rna_input_files = [os.path.join(folder, filename)
                           for filename in rna_files]
        func = lambda file, *_: int(run_get_stdout(f'wc -l < {file}', shell=True)) - 1
        self.run_function(func, folder, 
                          rna_input_files, repeat('', len(rna_input_files)))

    def run(self):
        """Calculate number of surviving pairs
        after each step"""
        # count in fastq folders
        if self.how == 'skip':
            logger.debug('Skipping read statistics calculation step')
            return
        for folder in ('rsites', 'dedup', 'trim'):
            self.count_in_fastqs(folder)
        # count BAM statistics
        if self.how == 'full':
            self.count_in_bams('hisat', 'mapped')
            self.count_in_bams('hisat', 'mapped2mism')
        # count in contacts
        self.count_contacts('contacts')
        # save result
        output_name = os.path.join(self.output_dir, f'{self.prefix}.tsv')
        result = pd.DataFrame.from_dict(self.result).sort_index()
        result.to_csv(output_name, sep='\t',
                      index=True, header=True)
