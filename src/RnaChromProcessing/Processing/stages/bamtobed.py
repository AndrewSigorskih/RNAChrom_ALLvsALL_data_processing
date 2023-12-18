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
        self.run_function(self._bam_to_bed, samples, output_samples)
        # return results
        return output_samples
    
    def _bam_to_bed(self,
                    inp_sample: SamplePair,
                    out_sample: SamplePair) -> int:
        """actual bam to bed conversion function"""
        exit_codes = []
        for infile, outfile in ((inp_sample.dna_file, out_sample.dna_file),
                                (inp_sample.rna_file, out_sample.rna_file)):
            cmd = (
                f"samtools view -Sh {infile} | grep -E 'NH:i:1$|^@' | samtools view -Sbh - |"
                f'bedtools bamtobed -cigar -i stdin > {outfile}'
            )
            exit_codes.append(run_command(cmd, shell=True))
        return exit_codes[0] or exit_codes[1]
