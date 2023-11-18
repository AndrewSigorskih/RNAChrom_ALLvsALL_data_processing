from pathlib import Path
from typing import List

from pydantic import PositiveInt

from .basicstage import BasicStage, SamplePair
from ...utils import run_command


class BamFilter(BasicStage):
    max_mismatch: PositiveInt = 2

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
    
    def run(self, samples: List[SamplePair]) -> None:
        """filter aligned reads: allow only reads that are aligned once
        with N or less mismatches"""
        # prepare filepaths
        dna_outputs = [
            self._stage_dir / sample.dna_file.name for sample in samples
        ]
        rna_outputs = [
            self._stage_dir / sample.rna_file.name for sample in samples
        ]
        # run function
        self.run_function(
            self._filter_bam,
            [sample.dna_file for sample in samples],
            [sample.rna_file for sample in samples],
            dna_outputs, rna_outputs
        )
        # save paths
        for sample, dna_out, rna_out in zip(samples, dna_outputs, rna_outputs):
            sample.set_files(dna_out, rna_out)
    
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
                f"grep -E 'XM:i:[0-{self.max_mismatch}]\s.*NH:i:1$|^@' | "
                f'samtools view -Sbh - > {outfile}'
            )
            exit_codes.append(run_command(cmd, shell=True))
        return exit_codes[0] or exit_codes[1]