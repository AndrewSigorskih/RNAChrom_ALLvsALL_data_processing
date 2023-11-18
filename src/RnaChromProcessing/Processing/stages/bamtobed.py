from pathlib import Path
from typing import List

from .basicstage import BasicStage, SamplePair
from ...utils import run_command


class BamToBed(BasicStage):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
    
    def run(self, samples: List[SamplePair]) -> List[SamplePair]:
        """bam to bed conversion"""
        # prepare filepaths
        output_samples = self._make_output_samples(samples, new_suff='bed')
        # run function
        self.run_function(
            self._bam_to_bed,
            [sample.dna_file for sample in samples],
            [sample.rna_file for sample in samples],
            [sample.dna_file for sample in output_samples],
            [sample.rna_file for sample in output_samples]
        )
        # return results
        return output_samples
    
    def _bam_to_bed(self,
                    dna_in_file: Path,
                    rna_in_file: Path,
                    dna_out_file: Path,
                    rna_out_file: Path) -> int:
        """actual bam to bed conversion function"""
        exit_codes = []
        for infile, outfile in ((dna_in_file, dna_out_file),
                                (rna_in_file, rna_out_file)):
            cmd = (
                f'samtools view -Sbh {infile} | '
                f'bedtools bamtobed -cigar -i stdin > {outfile}'
            )
            exit_codes.append(run_command(cmd, shell=True))
        return exit_codes[0] or exit_codes[1]
