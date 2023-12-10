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
        self.run_function(self._filter_bam, samples, output_samples)
        # return results
        return output_samples
    
    def _filter_bam(self,
                    inp_sample: SamplePair,
                    out_sample: SamplePair) -> int:
        """actual bam-filtering function"""
        exit_codes = []
        for infile, outfile in ((inp_sample.dna_file, out_sample.dna_file),
                                (inp_sample.rna_file, out_sample.rna_file)):
            cmd = (
                f'samtools view -Sh -F 4 {infile} | '
                f"grep -E 'XM:i:[0-{self.max_mismatch}]\s.*NH:i:1$|^@' | "
#                f"grep -E 'XM:i:[0-{self.max_mismatch}]|^@' | "
                f'samtools view -Sbh - > {outfile}'
            )
            exit_codes.append(run_command(cmd, shell=True))
        return exit_codes[0] or exit_codes[1]
