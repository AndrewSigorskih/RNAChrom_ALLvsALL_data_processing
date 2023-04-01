import concurrent.futures
import os
from typing import List

import pandas as pd

from .basicstage import suff_to_filter
from ...utils import find_in_list, run_command, run_get_stdout

class StatsCalc:
    def __init__(self, 
                 output_dir: str,
                 cpus: int,
                 dna_ids: List[str],
                 rna_ids: List[str]):
        self.output_dir: str = output_dir
        self.cpus: int = cpus
        self.dna_ids = dna_ids
        self.rna_ids = rna_ids
        self.result = {}
        os.makedirs('stats', exist_ok=True)

    def count_in_fastq_pair(self,
                            dna_file: str,
                            rna_file: str) -> int:
        """Counts matching IDs in pair of FASTQ files.
        Supported IDs formats:
            * @IDXXX @IDXXX\n
            * @FILE1.IDXXX @FILE2.IDXXX\n
            * @FILE.IDXXX @FILE.IDXXX"""
        dna_tmp_file = os.path.join('stats', rna_file.split('/')[-1] + '.tmp')
        rna_tmp_file = os.path.join('stats', dna_file.split('/')[-1] + '.tmp')

        for infile, outfile in ((rna_file, rna_tmp_file), (dna_file, dna_tmp_file)):
            cat = 'zcat' if infile.endswith('.gz') else 'cat'
            cmd = (
                f'{cat} {infile} | sed -n  "1~4p" | sed "s/@//g"  | '
                f'cut -d " " -f 1 | cut -d "." -f 2 > {outfile}'
            )
            run_command(cmd, shell=True)

        count_cmd =(
            "awk 'FNR==NR{array[$0]; next} ($0 in array) { count++ } END { print count } ' "
            f"{dna_tmp_file} {rna_tmp_file}"
        )
        result = run_get_stdout(count_cmd, shell=True)
        return int(result)
    
    def count_in_bam_pair(self,
                          dna_file: str,
                          rna_file: str,
                          mode='mapped') -> int:
        """Counts matching IDs in pair od BAM files.
        Supported IDs formats:
            * @IDXXX @IDXXX\n
            * @FILE1.IDXXX @FILE2.IDXXX\n
            * @FILE.IDXXX @FILE.IDXXX"""
        dna_tmp_file = os.path.join('stats', rna_file.split('/')[-1] + '.tmp')
        rna_tmp_file = os.path.join('stats', dna_file.split('/')[-1] + '.tmp')
        for infile, outfile in ((rna_file, rna_tmp_file), (dna_file, dna_tmp_file)):
            if mode == 'mapped':
                cmd = (
                    f"samtools view -F 4 {infile} | "
                    "awk -F \"\\t\" '{print \"@\"$1}' | "
                    "cut -d '.' -f2 | awk '{ a[$1]++ } END { for (b in a) { print b } }' "
                    f"> {outfile}"
                )
            elif mode == 'mapped2mism':
                cmd = (
                    f"samtools view -F 4 {infile} | grep -E 'XM:i:[0-2]\s.*' | "
                    "awk -F \"\\t\" '{print \"@\"$1}' | cut -d '.' -f2 | "
                    "awk '{ a[$1]++ } END { for (b in a) { print b } }' "
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
        with concurrent.futures.ThreadPoolExecutor(max_workers=self.cpus) as executor:
            future_to_id = {
                executor.submit(self.count_in_fastq_pair, dna_file, rna_file) : rna_id
                for dna_file, rna_file, rna_id in zip(dna_input_files, rna_input_files, self.rna_ids)
            }
            for future in concurrent.futures.as_completed(future_to_id):
                rna_id = future_to_id[future]
                count = future.result()
                self.result[folder][rna_id] = count
    
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
        with concurrent.futures.ThreadPoolExecutor(max_workers=self.cpus) as executor:
            future_to_id = {
                executor.submit(self.count_in_bam_pair, dna_file, rna_file, mode) : rna_id
                for dna_file, rna_file, rna_id in zip(dna_input_files, rna_input_files, self.rna_ids)
            }
            for future in concurrent.futures.as_completed(future_to_id):
                rna_id = future_to_id[future]
                count = future.result()
                self.result[mode][rna_id] = count
        
    def count_contacts(self, folder: str) -> None:
        self.result[folder] = {}
        filenames: List[str] = [x for x in os.listdir(folder)]
        rna_files = (find_in_list(id, filenames) for id in self.rna_ids)
        rna_input_files = (os.path.join(folder, filename)
                           for filename in rna_files)
        func = lambda file: int(run_get_stdout(f'wc -l < {file}', shell=True)) - 1
        with concurrent.futures.ThreadPoolExecutor(max_workers=self.cpus) as executor:
            future_to_id = {
                executor.submit(func, rna_file) : rna_id
                for rna_file, rna_id in zip(rna_input_files, self.rna_ids)
            }
            for future in concurrent.futures.as_completed(future_to_id):
                rna_id = future_to_id[future]
                count = future.result()
                self.result[folder][rna_id] = count

    def run(self):
        """Calculate number of surviving pairs
        after each step"""
        # count in fastq folders
        for folder in ('dedup', 'rsites', 'trim'):
            self.count_in_fastqs(folder)
        # count BAM statistics
        self.count_in_bams('hisat', 'mapped')
        self.count_in_bams('hisat', 'mapped2mism')
        # count in contacts
        self.count_contacts('contacts')
        # save result
        output_name = os.path.join(self.output_dir, 'stats.tsv')
        result = pd.DataFrame.from_dict(self.result).sort_index()
        result.to_csv(output_name, sep='\t',
                      index=True, header=True)
