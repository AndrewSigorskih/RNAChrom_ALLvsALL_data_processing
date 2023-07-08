from typing import List

from .basicstage import BasicStage
from ...utils import run_command

class AlignedToBed(BasicStage):
    def __init__(self, 
                 *args, **kwargs):
        super().__init__(*args, **kwargs)
    
    def run(self,
            dna_ids: List[str],
            rna_ids: List[str]):
        """bam to bed conversion"""
        func = self._aligned_to_bed
        self.run_function(func, dna_ids, rna_ids)

    def _aligned_to_bed(self,
                        dna_in_file: str,
                        rna_in_file: str,
                        dna_out_file: str,
                        rna_out_file: str):
        """actual bam to bed conversion function"""
        # change file extensions from NAME.bam to NAME.bed
        dna_out_file = dna_out_file.rsplit('.', 1)[0] + '.bed'
        rna_out_file = rna_out_file.rsplit('.', 1)[0] + '.bed'
        exit_codes = []
        for infile, outfile in ((dna_in_file, dna_out_file),
                                (rna_in_file, rna_out_file)):
            cmd = (
                f'samtools view -Sbh {infile} | '
                f'bedtools bamtobed -cigar -i stdin > {outfile}'
            )
            exit_code = run_command(cmd, shell=True)
            exit_codes.append(exit_code)
        return exit_codes[0] or exit_codes[1]