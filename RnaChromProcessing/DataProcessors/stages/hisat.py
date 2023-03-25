import shutil

from typing import Any, Dict, List

from .basicstage import BasicStage
from ...utils import check_file_exists, exit_with_error, run_command


class Hisat(BasicStage):
    def __init__(self, 
                 cfg: Dict[str, Any],
                 *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.tool: str = cfg.get('tool', 'hisat2')
        self.tool_path: str = cfg.get('tool_path', shutil.which(self.tool))
        self.genome_path: str = cfg.get('genome', None)
        self.known_splice: str = cfg.get('known_splice', None)
        self.hisat_threads: int = cfg.get('hisat_threads', 1)
        if (cpus := cfg.get('cpus', None)):
            self.cpus: int = cpus
        if not self.tool_path:
            exit_with_error(f'Cannot deduce path to {self.tool} executable!')
        if (not self.genome_path) or (not self.known_splice):
            exit_with_error('Genome or known splice were not provided!')
        check_file_exists(self.known_splice)

    def run(self,
            dna_ids: List[str],
            rna_ids: List[str]):
        """align reads using hisat2 tool"""
        if self.tool == 'hisat2':
            func = self._run_hisat
        else:
            exit_with_error(f'Unknown alignment tool: {self.tool}!')
        self.run_function(func, dna_ids, rna_ids)

    def _run_hisat(self,
                   dna_in_file: str,
                   rna_in_file: str,
                   dna_out_file: str,
                   rna_out_file: str):
        """run hisat2"""
        # change file extensions from NAME.* or NAME.*.gz to NAME.bam
        dna_out_file = dna_out_file.rstrip('.gz').rsplit('.', 1)[0]+ '.bam'
        rna_out_file = rna_out_file.rstrip('.gz').rsplit('.', 1)[0]+ '.bam'
        # create cmd and run
        dna_cmd = (
            f'{self.tool_path} -x {self.genome_path} -p {self.hisat_threads} '
            f'--no-spliced-alignment -k 100 --no-softclip -U {dna_in_file} | '
            f'samtools view -bSh > {dna_out_file}'
        )
        rna_cmd = (
            f'{self.tool_path} -x {self.genome_path} -p {self.hisat_threads} -k 100 '
            f'--no-softclip --known-splicesite-infile {self.known_splice} --dta-cufflinks '
            f'--novel-splicesite-outfile {rna_out_file}.novel_splice '
            f'-U {rna_in_file} | samtools view -bSh > {rna_out_file}'
        )
        return_code_1 = run_command(dna_cmd, shell=True)
        return_code_2 = run_command(rna_cmd, shell=True)
        return (return_code_1 or return_code_2)


