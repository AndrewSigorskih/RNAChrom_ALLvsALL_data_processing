from pathlib import Path
from typing import List

from pydantic import NonNegativeInt

from .basicstage import BasicStage, SamplePair
from ...utils import run_command


class BamFilter(BasicStage):
    max_mismatch: NonNegativeInt = 2

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
    
    def run(self, samples: List[SamplePair]) -> List[SamplePair]:
        """filter aligned reads: allow only reads that are aligned once
        with N or less mismatches"""
        # prepare filepaths
        output_samples = self._make_output_samples(samples)
        # run function
        self.run_function(
            self._filter_bam,
            [sample.dna_file for sample in samples],
            [sample.rna_file for sample in samples],
            [sample.dna_file for sample in output_samples],
            [sample.rna_file for sample in output_samples]
        )
        # return results
        return output_samples
    
    def _filter_bam(self,
                    dna_in_file: Path,
                    rna_in_file: Path,
                    dna_out_file: Path,
                    rna_out_file: Path) -> int:
        """actual bam-filtering function"""
        exit_codes = []
        for infile, outfile in ((dna_in_file, dna_out_file),
                                (rna_in_file, rna_out_file)):
            cmd = (
                f'samtools view -Sh -F 4 {infile} | '
#                f"grep -E 'XM:i:[0-{self.max_mismatch}]\s.*NH:i:1$|^@' | "
                f"grep -E 'XM:i:[0-{self.max_mismatch}]|^@' | "
                f'samtools view -Sbh - > {outfile}'
            )
            exit_codes.append(run_command(cmd, shell=True))
        return exit_codes[0] or exit_codes[1]
