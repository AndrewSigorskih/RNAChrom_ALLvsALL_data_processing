import os
import json
import argparse
import subprocess

parser = argparse.ArgumentParser(
    description=("Counts matching contacts in RNA and DNA files. "
                 "Resulting table will be stored in .tsv format in file "
                 "{BASE_DIR}/stats_{MODE}.lst (see args)"),
    epilog="requires: samtools>=1.15.1, grep, cut, wc, awk"
)

parser.add_argument("--config", required=True, help=("Path to config file in json format that contains "
                                                     "key 'base_dir' that specifies directory to run analysis in and "
                                                     "keys 'dna_ids' and 'rna_ids' corresponding to matching "
                                                     "ids of RNA and DNA files, e.g.:"
                                                     "dna_ids: [file1.dna, file2.dna, ...]"
                                                     "rna_ids: [file1.rna, file2.rna, ...]"
                                                   ))
parser.add_argument("--mode", required=False, help=("Specify filtration mode applied to aligned reads before counting. "
                                                    "'mapped': only omit unmapped reads i.e. samtools view -F 4, "
                                                    "'mapped2mism': only count reads that were mapped with 2 or less mismatches, "
                                                    "'unique': run for reads that were already filtered by unuqueness and number of mismatches"), 
                    choices=["mapped", "mapped2mism", "unique"], default="mapped")
parser.add_argument("--subdir", required=False, help=("Subdirectory of BASE_DIR to take files from."
                                                    "(Default: 'sam/')"), default="sam")
args = parser.parse_args()
with open(args.config, 'r') as f:
    dct = json.load(f)
    RNA_IDS = dct["rna_ids"]
    DNA_IDS = dct["dna_ids"]
    BASEDIR = dct["base_dir"]
    del dct
    
OUTFILENAME = os.path.join(BASEDIR, f"stats_{args.mode}.lst")
with open(OUTFILENAME, "w") as f:
    print("id1\tid2\ttot_reads_1\ttot_reads_2\tidentical", file=f)

    
def run_command(command: str):
    process = subprocess.Popen(command,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE,
                                shell=True,
                                executable="/bin/bash")
    (stdout, stderr) = process.communicate()
    if stderr:
        print(stderr.decode())
        return None
    return stdout.decode().strip()

def count_pairs(name1: str, name2: str):
    first = os.path.join(BASEDIR, f"tmp_{name1}.txt")
    second = os.path.join(BASEDIR, f"tmp_{name2}.txt")
    cnt_1 = run_command(f"wc -l < {first}")
    cnt_2 = run_command(f"wc -l < {second}")
    res = run_command(f"awk 'a[$0]++' {first} {second} | wc -l")
    run_command(f"echo -e {name1}\t{name2}\t{cnt_1}\t{cnt_2}\t{res} >> {OUTFILENAME}")
    
def count_aligned_pairs(name1, name2):
    for name in (name1, name2):
        infile = os.path.join(BASEDIR, args.subdir, f"{name}.bam")
        outfile = os.path.join(BASEDIR, f"tmp_{name}.txt")
        mode = args.mode
        if (mode == 'unique'):
            command = (f"samtools view {infile} | "
                       " awk -F \"\\t\" '{print \"@\"$1}' | "
                       f"cut -d'.' -f2 > {outfile}"
                      )
        elif (mode == 'mapped2mism'):
            command = (f"samtools view -F 4 {infile} | grep -E 'XM:i:[0-2]\s.*' | "
                       "awk -F \"\\t\" '{print \"@\"$1}' | cut -d'.' -f2 | "
                       "awk '{ a[$1]++ } END { for (b in a) { print b } }' " 
                       f" > {outfile}"
                      )
        elif (mode == 'mapped'):
            command = (f"samtools view -F 4 {infile} | "
                       " awk -F \"\\t\" '{print \"@\"$1}' | "
                       "cut -d'.' -f2 | awk '{ a[$1]++ } END { for (b in a) { print b } }' "
                       f"> {outfile}"
                      )
        else:
            print('Unknown mode, exiting!')
            exit()
        run_command(command)
    count_pairs(name1, name2)
    run_command(f"rm {os.path.join(BASEDIR, 'tmp_*.txt')}")

for rna_id, dna_id in zip(RNA_IDS, DNA_IDS):
    count_aligned_pairs(rna_id, dna_id)
