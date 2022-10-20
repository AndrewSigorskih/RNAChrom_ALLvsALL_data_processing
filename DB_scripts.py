import subprocess
import numpy as np
import pandas as pd
import itertools
import re

def run_command(command):
    process = subprocess.Popen(command,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE,
                                shell=True,
                                executable="/bin/bash")
    (stdout, stderr) = process.communicate()
    if stderr:
        print(stderr)
        return stderr
        #return stderr.decode()
    else:
        return stdout
        #return stdout.decode()
    
def run_trimmomatic(file1, file2, outfile1, outfile2,
                        trimmomatic_path="/mnt/scratch/rnachrom/tools/Trimmomatic-0.39/trimmomatic-0.39.jar",
                        window=5, qual_th=26, minlen=15):

    command = "java -jar {0} PE -phred33 {1} {2} {3} {4} {5} {6} SLIDINGWINDOW:{7}:{8} MINLEN:{9}".format(
            trimmomatic_path, file1, file2,
            outfile1, outfile1+".unpaired",
            outfile2, outfile2+".unpaired",
            int(window), int(qual_th), int(minlen))
    run_command(command)
    return None

def align(infile, outfile, mode,
         genome='/mnt/scratch/rnachrom/data/genomes/hg38/hg38',
         nthreads=8,
         novel_splice="",
         known_splice="/mnt/scratch/rnachrom/data/genes/human/splicesites.txt",
         bin_hisat2="/mnt/scratch/rnachrom/tools/hisat2/bin/hisat2",
         custom_command=None):
    
    if isinstance(custom_command, str):
        command=custom_command
    elif mode=="dna":
        command = r"""{} -p {} -x {} --no-spliced-alignment -k 100 --no-softclip """ \
                      r"""-U {} | samtools view -bSh > {}""".format(bin_hisat2, nthreads, genome, infile, outfile)
    elif mode=='rna':
        command = r"""{} -p {} -x {}  -k 100 --no-softclip --known-splicesite-infile {} --dta-cufflinks """ \
            r"""--novel-splicesite-outfile {} -U {} | samtools view -bSh > {}""".format(bin_hisat2, nthreads, genome, known_splice,
                                                                             novel_splice, infile, outfile)
    else:
        raise Exception("Mode unknown!")
    run_command(command)
    return None

def filter_aligned(infile, outfile):
    command = r"""samtools view -Sh -F 4 {} | grep -E 'XM:i:[0-2]\s.*NH:i:1$|^@' | samtools view -Sbh - > {}""".format(infile, outfile)
    run_command(command)
    return None

def aligned_to_bed(infile, outfile, mode='full'):
    if mode=='full':
        command = r"""samtools view -Sh -F 4 {} | grep -E 'XM:i:[0-2]\s.*NH:i:1$|^@' | samtools view -Sbh - | bedtools bamtobed -cigar -i stdin > {}""".format(infile, outfile)
    elif mode=='mapped2mism':
        command = r"""samtools view -Sh -F 4 {} | grep -E 'XM:i:[0-2]\s.*' | samtools view -Sbh - | bedtools bamtobed -cigar -i stdin > {}""".format(infile, outfile)
    elif mode == 'mapped':
        command = r"""samtools view -Sh -F 4 {} | samtools view -Sbh - | bedtools bamtobed -cigar -i stdin > {}""".format(infile, outfile)
    else:
        print('unknown mode, exiting!')
        return 9
    run_command(command)
    return None


def extract_contacts(infile_dna, infile_rna, outfile, chr_dict=False):
    dna = pd.read_csv(infile_dna, sep='\t', header=None)
    rna = pd.read_csv(infile_rna, sep='\t', header=None)
    dna[3] = dna[3].apply(lambda x: x.split('.')[1])
    rna[3] = rna[3].apply(lambda x: x.split('.')[1])
    dna.set_index([3], inplace=True)
    rna.set_index([3], inplace=True)
    
    ids = list(frozenset(dna.index).intersection(rna.index))
    
    dna = dna.loc[ids,]
    rna = rna.loc[ids,]
    
    res = pd.DataFrame(index=range(len(ids)),\
                       columns=['id', 'rna_chr', 'rna_bgn', 'rna_end', 'rna_strand', 'rna_cigar',\
                               'dna_chr', 'dna_bgn', 'dna_end', 'dna_strand', 'dna_cigar'])
    for i,item in enumerate(ids):
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
        
    if chr_dict:
        res['rna_chr'] = res['rna_chr'].apply(lambda x: chr_dict[x])
        res['dna_chr'] = res['dna_chr'].apply(lambda x: chr_dict[x])
    res.to_csv(outfile, sep='\t', index=False, header=True)


def extract_contacts_old(ids_rna, ids_dna, indir, outdir):
    for rnaf, dnaf in zip(ids_rna, ids_dna):
        dna = pd.read_csv('{}/{}.dna.bed'.format(indir, dnaf),sep='\t', header=None)
        dna[3] = dna[3].apply(lambda x: x.split('.')[1])
        dna.set_index([3], inplace=True)
        rna = pd.read_csv('{}/{}.rna.bed'.format(indir, rnaf),sep='\t', header=None)
        rna[3] = rna[3].apply(lambda x: x.split('.')[1])
        rna.set_index([3], inplace=True)
        d = {'read_id':[], 'rna_chr':[], 'rna_start':[], 'rna_end':[], 'rna_strand':[],\
         'rna_cigar':[], 'rna_nlen':[], 'dna_chr':[], 'dna_start':[], 'dna_end':[],\
         'dna_strand':[], 'dna_cigar':[]}
        for item in dna.index:
            if item in rna.index:
                d['read_id'].append(item)
                d['rna_chr'].append(rna.at[item,0])
                d['rna_start'].append(rna.at[item,1])
                d['rna_end'].append(rna.at[item,2])
                d['rna_strand'].append(rna.at[item,5])
                cigar = rna.at[item,6]
                nlen = 0
                if ('N') in cigar:
                    nlen = re.split(r'\D+',cigar.split('N')[0])[-1]
                    nlen = int(nlen)
                d['rna_cigar'].append(cigar)
                d['rna_nlen'].append(nlen)
                d['dna_chr'].append(dna.at[item,0])
                d['dna_start'].append(dna.at[item,1])
                d['dna_end'].append(dna.at[item,2])
                d['dna_strand'].append(dna.at[item,5])
                d['dna_cigar'].append(dna.at[item,6])
        data = pd.DataFrame(d)
        data.to_csv('{}/{}.bed'.format(outdir, rnaf), sep='\t', index=False, header=True)
        print('{} done!'.format(rnaf))