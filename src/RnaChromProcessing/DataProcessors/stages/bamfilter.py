from typing import List

from .basicstage import BasicStage
from ...utils import run_command

class BamFilter(BasicStage):
    def __init__(self, 
                 *args, **kwargs):
        super().__init__(*args, **kwargs)
    
    def run(self,
            dna_ids: List[str],
            rna_ids: List[str]):
        """filter aligned reads: allow only reads that are aligned once
        with 2 or less mismatches"""
        func = self._filter_bam
        self.run_function(func, dna_ids, rna_ids)
    
    def _filter_bam(self,
                    dna_in_file: str,
                    rna_in_file: str,
                    dna_out_file: str,
                    rna_out_file: str):
        """actual bam-filtering function"""
        exit_codes = []
        for infile, outfile in ((dna_in_file, dna_out_file),
                                (rna_in_file, rna_out_file)):
            cmd = (
                f'samtools view -Sh -F 4 {infile} | '
                "grep -E 'XM:i:[0-2]\s.*NH:i:1$|^@' | "
                f'samtools view -Sbh - > {outfile}'
            )
            exit_code = run_command(cmd, shell=True)
            exit_codes.append(exit_code)
        return exit_codes[0] or exit_codes[1]
