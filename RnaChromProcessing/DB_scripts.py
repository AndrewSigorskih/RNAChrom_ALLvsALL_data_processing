import os
import shutil
import subprocess
from tempfile import NamedTemporaryFile
from typing import Optional

import numpy as np
import pandas as pd

def _run_command(command: str) -> str:
    print(command)
    process = subprocess.Popen(command,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE,
                                shell=True,
                                executable="/bin/bash")
    (stdout, stderr) = process.communicate()
    return stdout.decode(), stderr.decode()

    
def _run_check_command(command: str):
    (stdout, stderr) = _run_command(command)
    if stderr:
        print(stderr)
    return None

class DB_processor:
    def __init__(self, cfg_file: str):
        import json
        with open(cfg_file, "r") as f:
            dct = json.load(f)
        self.base_dir = dct.get("base_dir", os.getcwd())
        self.fastuniq_pth = dct.get("fastuniq_pth", shutil.which('fastuniq'))
        self.trimmomatic_pth = dct.get("trimmomatic_pth", shutil.which('trimmomatic'))
        self.hisat2_pth = dct.get("hisat2_pth", "hisat2")
        self.genome = dct.get("genome_pth", "/mnt/scratch/rnachrom/data/genomes/hg38/hg38")
        self.known_splice = dct.get("known_splice","/mnt/scratch/rnachrom/data/genes/human/splicesites.txt")
        self.chr_dct = dct.get("chr_dct", None)
        self.nthreads = dct.get("nthreads", 8)
        self.rna_ids = dct.get('rna_ids', None)
        self.dna_ids = dct.get('dna_ids', None)
        if (not self.rna_ids) or (not self.dna_ids):
            raise ValueError("File ids not specified!")
        # todo: manage dirs
        # todo: add checks for non-existing paths to all functions
    
    def run_fastuiniq(self,
                      indir: Optional[str] = None,
                      outdir: Optional[str] = None):
        """delete non-unique reads"""
        if not indir: 
            indir = os.path.join(self.base_dir, "raw")
        if not outdir:
            outdir = os.path.join(self.base_dir, "fastuniq")
        for dna, rna in zip(self.dna_ids, self.rna_ids):
            infile1 = os.path.join(indir, f'{dna}.fastq')
            infile2 = os.path.join(indir, f'{rna}.fastq')
            outfile1 = os.path.join(outdir, f'{dna}.fastq')
            outfile2 = os.path.join(outdir, f'{rna}.fastq')
            with NamedTemporaryFile(mode='w', delete=False, dir='.') as f:
                tmpfilename = f.name
                print(infile1, infile2, sep='\n', file=f)
            command = (f'{self.fastuniq_pth} -i {tmpfilename} '
                       f'-o {outfile1} -p {outfile2}')
            _run_check_command(command)
            os.remove(tmpfilename)
                   
    def run_trimmomatic(self,
                        window: int = 5,
                        qual_th: int = 26,
                        minlen: int = 15,
                        indir: Optional[str] = None,
                        outdir: Optional[str] = None):
        """trim reads by quality"""
        if not indir: 
            indir = os.path.join(self.base_dir, "fastuniq")
        if not outdir:
            outdir = os.path.join(self.base_dir, "trimmed")
        for dna, rna in zip(self.dna_ids, self.rna_ids):
            command = (f"{self.trimmomatic_pth} PE -phred33 "
                       f"{os.path.join(indir, f'{rna}.fastq')} {os.path.join(indir, f'{dna}.fastq')} "
                       f"{os.path.join(outdir, f'{rna}.fastq')} {os.path.join(outdir, f'{rna}.unpaired')} "
                       f"{os.path.join(outdir, f'{dna}.fastq')} {os.path.join(outdir, f'{dna}.unpaired')} "
                       f"SLIDINGWINDOW:{window}:{qual_th} MINLEN:{minlen}"
            )
            _run_check_command(command)
        return None
    
    def align(self):
        """align reads to genome"""
        indir = os.path.join(self.base_dir, "trimmed")
        outdir = os.path.join(self.base_dir, "sam")
        for filename in self.dna_ids:
            infilename = os.path.join(indir, f'{filename}.fastq')
            outfilename = os.path.join(outdir, f'{filename}.bam')
            command = (f"{self.hisat2_pth} -p {self.nthreads} -x {self.genome} "
                       f"--no-spliced-alignment -k 100 --no-softclip -U "
                       f"{infilename} | "
                       f"samtools view -bSh > {outfilename}"          
            )
            _run_check_command(command)
            
        for filename in self.rna_ids:
            infilename = os.path.join(indir, f'{filename}.fastq')
            outfilename = os.path.join(outdir, f'{filename}.bam')
            command = (f"{self.hisat2_pth} -p {self.nthreads} -x {self.genome} "
                       f"-k 100 --no-softclip --known-splicesite-infile {self.known_splice} "
                       f"--dta-cufflinks --novel-splicesite-outfile {outfilename}.novel_splice "
                       f"-U {infilename} | samtools view -bSh > "
                       f"{outfilename}"
            )
            _run_check_command(command)
        return None
    
    def filter_aligned(self):
        """filter aligned reads: allow only readt that are aligned once
        with 2 or less mismatches"""
        for filename in self.rna_ids + self.dna_ids:
            infilename = os.path.join(self.base_dir, "sam", f"{filename}.bam")
            outfilename = os.path.join(self.base_dir, "bam", f"{filename}.bam")
            command = (f"samtools view -Sh -F 4 {infilename} | "
                       f"grep -E 'XM:i:[0-2]\s.*NH:i:1$|^@' | "
                       f"samtools view -Sbh - > {outfilename}"
                      )
            _run_check_command(command)
        return None
        
    
    def aligned_to_bed(self,
                       indir: str = "bam",
                       outdir: str = "bed",
                       mode: str = "full"):
        """bam to bed conversion"""
        for filename in self.rna_ids + self.dna_ids:
            infile = os.path.join(self.base_dir, indir, f"{filename}.bam")
            outfile = os.path.join(self.base_dir, outdir, f"{filename}.bed")
            if mode=='full':
                command = (f"samtools view -Sh -F 4 {infile} | "
                           f"grep -E 'XM:i:[0-2]\s.*NH:i:1$|^@' | samtools view -Sbh - | "
                           f"bedtools bamtobed -cigar -i stdin > {outfile}"
                )
            elif mode == "mapped2mism":
                command = (f"samtools view -Sh -F 4 {infile} | "
                           f"grep -E 'XM:i:[0-2]\s.*' | samtools view -Sbh - | "
                           f"bedtools bamtobed -cigar -i stdin > {outfile}"
                )
            elif mode == "mapped":
                command = (f"samtools view -Sh -F 4 {infile} | "
                           f"samtools view -Sbh - | "
                           f"bedtools bamtobed -cigar -i stdin > {outfile}"
                )
            else:
                raise ValueError(f"Unknown mode: {mode}")
            _run_check_command(command)
        return None
    
    def _extract_contacts(self, 
                          infile_dna: str,
                          infile_rna: str,
                          outfile: str):
        """take 2 bed files, produce 1 combined contacts file"""
        dna = pd.read_csv(infile_dna, sep='\t', header=None)
        rna = pd.read_csv(infile_rna, sep='\t', header=None)
        dna[3] = dna[3].apply(lambda x: x.split('.')[1])
        rna[3] = rna[3].apply(lambda x: x.split('.')[1])
        dna.set_index([3], inplace=True)
        rna.set_index([3], inplace=True)
        #ids = list(frozenset(dna.index).intersection(rna.index))
        ids = list(set(dna.index) & set(rna.index))                  
        dna = dna.loc[ids,]
        rna = rna.loc[ids,]
        res = pd.DataFrame(index=range(len(ids)),\
                           columns=['id', 'rna_chr', 'rna_bgn', 'rna_end', 'rna_strand', 'rna_cigar',\
                                    'dna_chr', 'dna_bgn', 'dna_end', 'dna_strand', 'dna_cigar'])
        for i, item in enumerate(ids):
            res.at[i,'id'] = item
            res.at[i,'rna_chr'] = rna.at[item,0]
            res.at[i,'rna_bgn'] = rna.at[item,1]
            res.at[i,'rna_end'] = rna.at[item,2]
            res.at[i,'rna_strand'] = rna.at[item,5]
            res.at[i,'rna_cigar'] = rna.at[item,6]
            res.at[i,'dna_chr'] = dna.at[item,0]
            res.at[i,'dna_bgn'] = dna.at[item,1]
            res.at[i,'dna_end'] = dna.at[item,2]
            res.at[i,'dna_strand'] = dna.at[item,5]
            res.at[i,'dna_cigar'] = dna.at[item,6]
        
        if self.chr_dct:
            res['rna_chr'] = res['rna_chr'].apply(lambda x: self.chr_dict[x])
            res['dna_chr'] = res['dna_chr'].apply(lambda x: self.chr_dict[x])
        res.to_csv(outfile, sep='\t', index=False, header=True)
                           
    def make_contacts(self):
        """make contacts files for all pairs of ids"""
        for dna, rna in zip(self.dna_ids, self.rna_ids):
            self._extract_contacts(os.path.join(self.base_dir, "bed", f"{dna}.bed"),
                                   os.path.join(self.base_dir, "bed", f"{rna}.bed"),
                                   os.path.join(self.base_dir, "contacts", f"{rna}.tab")
            )
        return None
