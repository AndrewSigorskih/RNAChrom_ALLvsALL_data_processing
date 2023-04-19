# RNAChrom_ALLvsALL_data_processing

Pachage for processing data of several "ALL-vs-ALL" RNA-Chromatin interactions capturing experiments, e.g. RADICL, GRID, iMARGI, CHAR, etc.

## Table of Contents
1. [Overwiew](#overview)
2. [Installation](#installation)
    1. [Required software](#software1)
    2. [Setting up virtual environment](#venv)
    3. [Actual installation](#install1)
3. [Usage](#usage)
4. [X-RNA inference](#xrna)

<a name="overview"></a>
##  Overview
This package is desined to unify the data processing of different RNA-Chromatin interactions capturing experiments.
The main routine performs several consecutive steps of data processing as follows:

1. Dedupication step. PCR duplicates are removed from read pairs pool.
2. Restriction sites handling. Filtering reads by starting bases  and/or finishing the fragmented restriction site. Several "procedure persets" are available for this step based on previous papers, hovewer users can also use their own scripts by passing those as innputs (see "Config contents" in the "[Usage](#usage)" section for more information).
3. Trimming. Trimming of the input fastq files by quality. The default behaviour for this step is using Trimmomatic with parameters **SLIDINGWINDOW:5:26 MINLEN:15**.
4. Mapping RNA and DNA reads to reference genome separately using the Hisat2 software.
5. Filtering resulting bam files: only reads that were uniquely mapped with 2 or fewer mismatches will pass this stage.
6. Converting filtered bam files to bed and joining them into a resulting contacts table.
7. Overall statistics calculation: how many pairs of reads survived each step.

Each step from 1 to 3 can be skipped. The package interface provides opportunity to fine-tune the behaviour for almost each step, but for running defaults the user's input can be minimalistic.


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

## Usage

```
rnachromprocessing [-h] -c CONFIG [-mode {contacts,XRNA}] [-v]
```

Supported arguments:
```
-h, --help            show help message and exit
-c CONFIG, --config CONFIG
                    Configuration file (required).
-mode {contacts,XRNA}, --run_mode {contacts,XRNA}
                    Action to be performed during run:
                    * contacts: build contacts data from provided raw FASTQ files (default)
                    * XRNA: attempt to infer X-RNAs from provided contacts and BAM files.
-v, --verbose         Add logging verbosity.
```

All inputs and parameters are passed using the </i>config</i> argument -- a path to existing file in a <b>.json</b> format.

Config should contain several required fields and may contain several additional fields for fine-tuning. Complete config specification and examples are listed below.

### Config contents

"Global" options:
* "rna_ids" : non-empty array of strings. Must contain names of files with RNA parts of reads, without directory names or file extensions. For example, for the file named <i>SRR9201799.fastq.gz</i> id is SRR9201799.
* "dna_ids" : non-empty array of strings.  Must contain names of files with DNA parts of reads, synchronized with rna_ids array.
* "base_dir" : string. Optional. Directory to run analysis in (one or several temporal directories will be created in this directory). Path should exist if provided. If not provided, current working directory will be taken.
* "input_dir" : string. Required. Path to directory with input DNA and RNA files.
* "output_dir" : string. Required. Path to directory to store results in. Will be created if doesnt exist.
* "keep" : non-empty array of strings. Should contain names of steps that will have their results saved into "output_dir".
    Supported values (can be listed in any order): ["dedup", "rsites", "trim", "hisat", "bam", "bed", "contacts"].
* "cpus" : int > 0. Number of tasks to run simultaneously on each step. Default is 1.

Most of the steps require several options, that are grouped into corresponding sub-configs:

* "dedup" : sub-config for deduplication step. Supported values:
    * "tool" : what deduplication tool to use. Default is "fastuniq". Supported options:
        * "fastuniq" : use [fastuniq](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0052249) tool.
        * "skip" : do not perform deduplication step.
    * "tool_path" : Optional. Path to tool executable. Use if your tool is not available via PATH. If this option is not set, package will try to search for executable location using the `shutil.which()` function with the name of tool provided in previous option.
    * "cpus": int > 0. Optional. Number of tasks to run simultaneously for deduplication step. If not set, the "global" value will be used.

* "rsites" : sub-config for restriction-site handling step. Supported options:
    * "type" : what action to perform on read pairs during this step. Default is "skip". Supported values:
        * "imargi" : Save read pair if DNA read starts with CT or NT and remove first 2 bases from RNA reads. Reads in files should be synchronized.
        * "grid" : Save read pair if DNA read ends with AG. Add CT to the end of selected DNA reads (quality is copied from AG bases). Reads in files should be synchronized.
        * "custom" : Use custom script provided by user.
        * "skip" : do not perform this step.
    * "tool_path" : Optional. Will be only considered if "custom" type of action is chosen for previous option. Path to user's script that will manage all operations on restriction sites for read pairs. Script should be executable and take exactly 4 positional arguments: `dna-input-file`, `rna-input-file`, `dna-output-file`, `rna-output-file`.
    * "cpus": int > 0. Optional. Number of tasks to run simultaneously for this step. If not set, the "global" value will be used.

* "trim" : sub-config for trimming step. Supported options:
    * "tool" : what trimming tool to use. Default is "trimmomatic". Supported options:
        * "trimmomatic"
        * "skip" : do not perform trimming step.
    * "params" : Optional. Tool-specific parameters in form of a sub-config. For "Trimmomatic" tool supported parameers are "window", "qual_th" and "minlen" (int > 0 each).
    * "cpus": int > 0. Optional. Number of tasks to run simultaneously for this step. If not set, the "global" value will be used.

* "hisat" : sub-config for alignment step. Supported options:
    * "genome" : **Required**. The basename of the index for the reference genome (will be passed to hisat's `-x` option; see more in [reference manual](http://daehwankimlab.github.io/hisat2/manual/#:~:text=number%3E%7D%20%5B%2DS%20%3Chit%3E%5D-,Main%20arguments,-%2Dx%20%3Chisat2%2Didx))
    * "known_splice" : **Required**. Path to known-splice file (will be passed to hisat's `--known-splicesite-infile` option).
    * "tool_path" : Optional. Path to hisat2 main executable if not available via PATH.
    * "cpus": int > 0. Optional. Number of tasks to run simultaneously for this step. If not set, the "global" value will be used.
    * "hisat_threads" : int > 0. How many threads hisat will use **for every** simultaneously-run task. (Value to pass to hisat's `-p` option). Default 1.

### Minimal config example

```
{
    "rna_ids": ["SRR9201799", "SRR9201801"],
    "dna_ids": ["SRR9201800", "SRR9201802"],
    "base_dir" : ".",
    "input_dir": "../data/radicl/raw",
    "output_dir": "results",
    "cpus": 2,
    "keep" : ["contacts", "bam", "trim", "dedup"],
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
    "keep" : ["contacts", "bam"],
    "dedup": {
        "tool": "fastuniq",
        "tool_path": "/path/to/fastuniq",
        "params": [],
        "cpus": None
    },
    "rsites": {
        "cpus": None,
        "type": "skip"
    },
    "trim": {
        "tool" : "trimmomatic",
        "tool_path" : "path/to/trimmomatic/jar",
        "cpus": None,
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
        "cpus": None
    }
}
```
<a name="xrna"></a>
## X-RNA inference 

Work in progress: this functional is not properly implemented yet!
