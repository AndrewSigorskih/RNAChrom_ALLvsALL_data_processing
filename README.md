# RNAChrom_ALLvsALL_data_processing

Package for processing data of several "ALL-vs-ALL" RNA-Chromatin interactions capturing experiments, e.g. RADICL, GRID, iMARGI, CHAR, etc.

## Table of Contents
1. [Overview](#overview)
2. [Installation](#installation)
    1. [Required software](#software1)
    2. [Setting up virtual environment](#venv)
    3. [Actual installation](#install1)
    4. [Data preparation](#genomeprep)
3. [rnachromprocessing](#rnachromprocessing)
    1. [Usage](#processingusage)
    2. [Config contents](#processingconfig)
    3. [Minimal config example](#minimalconfig)
    4. [Advanced config example](#advancedconfig)
    5. [Outputs](#processingoutputs)
4. [detect-strand](#detect-strand)
    1. [Usage](#strandusage)
    2. [Config contents](#strandconfig)
5. [X-RNA inference](#xrna)

<a name="overview"></a>
##  Overview
This package is desined to unify the data processing of different RNA-Chromatin interactions capturing experiments.

After installation it provides access to several command-line tools:
* rnachromprocessing -- pipeline for uniform processing of "All-versus-All" RNA-Chromatin interactions data. 
* detect-strand -- program that checks whether orientation of RNA parts of contacts was inverted during sequencing. 
* infer-xrna -- pipeline that searches for most promising not-annotated transcripts in the experimental data.

<a name="software1"></a>
## Installation
### Required software

Basic dependencies:

* Python (v 3.9 or higher)
* Common unix command-line text-processing utilities: grep, sed, awk, cut, uniq etc.

This package also requires several bioinformatic tools to run.

* samtools ver 1.11 or higher (should be available via PATH variable)
* bedtools (should be available via PATH variable)
* hisat2
* any of the supported trimming tools (trimmomatic, ...)
* any of the supported deduplication tools (fastuniq, ...)

<a name="venv"></a>
### Setting up virtual environment (recommended)
It is highly recommended to install package using in a virtual environment. This can be achieved via [python's venv](https://docs.python.org/3/library/venv.html) or [conda envs](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html).<br> Moreover, all of the required bioinformatic tools can be installed via conda as well (but mind the tricky [conda samtools installation](https://www.biostars.org/p/455593/)!) 


<a name="install1"></a>
### Clone and install package
Clone repository from github and install:

```
git clone https://github.com/AndrewSigorskih/RNAChrom_ALLvsALL_data_processing.git
cd RNAChrom_ALLvsALL_data_processing
pip install .
```

<a name="genomeprep"></a>
### Data preparation.
The *rnachromprocessing* pipeline expects reads in paired files to be synchronized, i.e. for every index i read number i from rna.fastq file and read number i from dna.fastq file have same id. If you download data from SRA archive using fasterq-dump make sure to use `--split-3` flag.

Chromosome names in genome and gene annotation should be the same. We also recommend to remove mitochondrial and "Unplaced" sequences from genome fasta. This can be easily achieved using tools like [seqkit](https://bioinf.shenwei.me/seqkit/): for example, command `seqkit grep -vrp "^chrUn" file.fa > clean.fa` will  remove all unplaced chromosome fragments from genome file.

Hisat2 requires genome to be indexed. In order to build the index, use command `hisat2-build -p 16 genome.fa prefix`, where **prefix** is the prefix name for genome index files. In order to properly map the RNA parts you will also need the splicecite file, that can be obtained from genome annotation file using hisat's script: `hisat2_extract_splice_sites.py genome.gtf > genome.ss`. See [hisat2 manual](http://daehwankimlab.github.io/hisat2/howto/#building-indexes) for more information.

## rnachromprocessing

The main routine performs several consecutive steps of data processing as follows:

1. Restriction sites handling. Filtering reads by starting bases and/or finishing the fragmented restriction site. Several "procedure presets" are available for this step based on previous papers, however users can also use their own scripts by passing those as innputs (see "Config contents" section for more information).
2. Dedupication step. PCR duplicates are identified and removed.
3. Trimming. Trimming of the input fastq files by quality. The default behaviour for this step is using Trimmomatic with parameters **SLIDINGWINDOW:5:26 MINLEN:15**.
4. Mapping RNA and DNA reads to reference genome separately using the Hisat2 software.
5. Filtering resulting bam files: only reads that were uniquely mapped with 2 or fewer mismatches will pass this stage.
6. Converting filtered bam files to bed and joining them into a resulting contacts table.
7. Overall statistics calculation: how many pairs of reads survived each step.

Each step from 1 to 3 can be skipped. The pipeine interface provides opportunity to fine-tune the behaviour for almost each step, but for running defaults the user's input can be minimalistic.

<a name="processingusage"></a>
### Usage

```
rnachromprocessing [-h] -c CONFIG [-v]
```

Supported arguments:
```
-h, --help          show help message and exit
-c CONFIG, --config CONFIG
                    Configuration file in json format (required).
-v, --verbose       Log debug information about each step to stdout.
```

All inputs and parameters are passed using the </i>config</i> argument -- a path to existing file in a <b>.json</b> format.

Config should contain several required fields and may contain several additional fields for fine-tuning. Complete config specification and examples are listed below.

**Config examples can be found in `configs/rnachromprocessing` directory.**

<a name="processingconfig"></a>
### Config contents

"Global" options:
* "rna_ids" : non-empty array of strings. Must contain names of files with RNA parts of reads, without directory names or file extensions. For example, for the file named <i>SRR9201799.fastq.gz</i> id is SRR9201799.
* "dna_ids" : non-empty array of strings.  Must contain names of files with DNA parts of reads, synchronized with rna_ids array.
* "base_dir" : string. Optional. Directory to run analysis in (one or several temporal directories will be created in this directory). Path should exist if provided. If not provided, current working directory will be taken.
* "input_dir" : string. Required. Path to directory with input DNA and RNA files.
* "output_dir" : string. Required. Path to directory to store results in. Will be created if doesnt exist.
* "keep" : non-empty array of strings. Should contain names of steps that will have their results saved into "output_dir".
    Supported values (can be listed in any order): ["rsites", "dedup", "trim", "hisat", "bam", "bed", "contacts"]. By default only "bam" and "contacts" steps are saved.
* "cpus" : int > 0. Number of tasks to run simultaneously on each step. Default is 1.

Most of the steps require several options, that are grouped into corresponding sub-configs:

* "rsites" : sub-config for restriction-site handling step. Supported options:
    * "type" : what action to perform on read pairs during this step. Default is "skip". Supported values:
        * "imargi" : Save read pair if DNA read starts with CT or NT and remove first 2 bases from RNA reads. Reads in files should be synchronized.
        * "grid" : Save read pair if DNA read ends with **rsite_bgn** bases (default AG). Add **rsite_end** bases (default CT) to the end of selected DNA reads (quality is copied from **rsite_bgn** bases). Reads in files should be synchronized.
        * "custom" : Use custom script provided by user. 
        * "skip" : do not perform this step.
    * "rsite_bgn": beginning of restriction site that will be used in "grid"-like processing. Default "AG".
    * "rsite_end": Missing end of restriction site that will be used in "grid"-like processing. Default "CT".
    * "tool_path" : Optional. Will be only considered if "custom" type of action is chosen for "type" option. Path to user's script that will manage all operations on restriction sites for read pairs. Script should be executable and take exactly 4 positional arguments: `dna-input-file`, `rna-input-file`, `dna-output-file`, `rna-output-file`.
    * "cpus": int > 0. Optional. Number of tasks to run simultaneously for this step. If not set, the "global" value will be used.

* "dedup" : sub-config for deduplication step. Supported values:
    * "tool" : what deduplication tool to use. Default is "fastuniq". Supported options:
        * "fastuniq" : use [fastuniq](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0052249) tool.
        * "skip" : do not perform deduplication step.
    * "tool_path" : Optional. Path to tool executable. Use if your tool is not available via PATH. If this option is not set, package will try to search for executable location using the `shutil.which()` function with the name of tool provided in previous option.
    * "cpus": int > 0. Optional. Number of tasks to run simultaneously for deduplication step. If not set, the "global" value will be used.

* "trim" : sub-config for trimming step. Supported options:
    * "tool" : what trimming tool to use. Default is "trimmomatic". Supported options:
        * "trimmomatic"
        * "skip" : do not perform trimming step.
    * "params" : Optional. Tool-specific parameters in form of a sub-config. For "Trimmomatic" tool supported parameters are "window", "qual_th" and "minlen" (int > 0 each).
    * "tool_path" : Optional. Path to tool executable. Use if your tool is not available via PATH. If this option is not set, package will try to search for executable location using the `shutil.which()` function with the name of tool provided in "type" option.
    * "cpus": int > 0. Optional. Number of tasks to run simultaneously for this step. If not set, the "global" value will be used.

* "hisat" : sub-config for alignment step. Supported options:
    * "genome" : **Required**. Path ending with prefix of genome index files (will be passed to hisat's `-x` option; see more in [reference manual](http://daehwankimlab.github.io/hisat2/manual/#:~:text=number%3E%7D%20%5B%2DS%20%3Chit%3E%5D-,Main%20arguments,-%2Dx%20%3Chisat2%2Didx))
    * "known_splice" : **Required**. Path to known-splice file (will be passed to hisat's `--known-splicesite-infile` option).
    * "tool_path" : Optional. Path to hisat2 main executable if not available via PATH.
    * "cpus": int > 0. Optional. Number of tasks to run simultaneously for this step. If not set, the "global" value will be used.
    * "hisat_threads" : int > 0. How many threads hisat will use **for every** simultaneously-run task. (Value to pass to hisat's `-p` option). Default 1.

<a name="minimalconfig"></a>
### Minimal config example

```
{
    "rna_ids": ["SRR9201799", "SRR9201801"],
    "dna_ids": ["SRR9201800", "SRR9201802"],
    "base_dir" : ".",
    "input_dir": "../data/radicl/raw",
    "output_dir": "results",
    "cpus": 2,
    "dedup": {
        "tool": "fastuniq",
        "tool_path": "tools/fastuniq"
    },
    "rsites": {
        "type": "skip"
    },
    "trim": {
        "tool" : "trimmomatic"
    },
    "hisat": {
        "genome": "data/mm10/mm10",
        "known_splice": "data/genes/gencode.vM32.ss"
    }
}
```
<a name="advancedconfig"></a>
### Advanced config example

```
{
    "rna_ids": [],
    "dna_ids": [],
    "schema": "radicl",
    "base_dir" : ".",
    "input_dir": "./fastq",
    "output_dir": "results",
    "cpus": 5,
    "keep" : ["contacts", "bam", "rsites"],
    "dedup": {
        "tool": "fastuniq",
        "tool_path": "/path/to/fastuniq",
        "params": [],
        "cpus": 2
    },
    "rsites": {
        "type": "imargi"
    },
    "trim": {
        "tool" : "trimmomatic",
        "tool_path" : "path/to/trimmomatic/jar",
        "params": {
            "window": 5,
            "qual_th": 26,
            "minlen": 15
        }
    },
    "hisat": {
        "tool_path" : "path/to/hisat",
        "genome": "PATH/TO/GENOMEDIR",
        "known_splice": "PATH/TO/KNOWN/SPLICE",
    }
}
```

<a name="processingoutputs"></a>
### Outputs

Tables in tsv format. Columns:
* 'id': read pair unique ID.
* 'rna_chr', 'rna_bgn', 'rna_end', 'rna_strand', 'rna_cigar': RNA read information from corresponfing bed file.
* 'dna_chr', 'dna_bgn', 'dna_end', 'dna_strand', 'dna_cigar': DNA read information from corresponfing bed file.


## detect-strand

This program takes contacts files from previous step as input and checks whether orientation of RNA parts of contacts was inverted or lost during sequencing.

The following procedure is performed:
* User provides gene annotation file in GTF format and file with a list of selected genes. This file can contain gene names or gene ids from 9th field of gene annotation file.
* For each contacts file, a subset of contacts RNA parts aligning with selected genes is selected.
* For each gene, the correponding contacts pile is divided into two parts: the ones with their strand matching the gene strand and the ones from the opposite strand. Pairs of values for each file for each gene are stored in the `raw_counts.tsv` file.
* For each file, genes "vote" for strand: numbers of "wins" for gene strand and opposite strand are obtained. Cases when both strands had exactly the same number of contacts or when no contacts at all are ommited.
* Resulting wins and losses for each file are visualized as barplots and saved in .png and .svg formats.
* Resulting wins and losses for each file are saved in the `wins.tsv` file with following columns: 'same' holds the number of wins, 'anti' holds the number of losses and 'strand' holds the decision: "SAME" if the genes' strands won overwhelmingly, "ANTI" if the opposite strands won overwhelmingly and "UNKNOWN" otherwise.

<a name="strandusage"></a>
### Usage

```
detect-strand [-h] -c CONFIG [-v]
```
Supported arguments:
```
-h, --help          show help message and exit
-c CONFIG, --config CONFIG
                    Configuration file in json format (required).
-v, --verbose       Log debug information about each step to stdout.
```

**Examples of configs and genes lists for human and mouse datasets can be found in `configs/strand` directory.**

<a name="strandconfig"></a>
### Config contents

* "input_dir": path to input directory with contacts tables.
* "output_dir": path to directory to store outputs in. Will be created if not exist.
* "gtf_annotation": path to gene annotation file in GTF format.
* "genes_list": path to a file with selected genes names (or ids), one per line
* "prefix": prefix name for all output files. Default  is "strand".
* "exp_groups": sub-config with lists of input file ids divided by groups (experiment types, for example). See config examples for more information.

<a name="xrna"></a>
## infer-xrna

Work in progress: this functional is not properly implemented yet!
