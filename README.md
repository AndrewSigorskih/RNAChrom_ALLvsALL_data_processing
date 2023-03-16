# RNAChrom_ALLvsALL_data_processing

Pachage for processing data of several "ALL-vs-ALL" RNA-Chromatin interactions capturing experiments, e.g. RADICL, GRID, iMARGI, CHAR, etc.

## Table of Contents
1. [Installation](#installation)
    1. [Required software](#software1)
    2. [Setting up virtual environment](#venv)
    3. [Actual installation](#install1)
2. [Usage](#usage)
3. [X-RNA inference](#xrna)


## Installation
<a name="software1"></a>
### Required software

This package requires several bioinformatic tools to run.

* samtools ver 1.11 or higher (should be available via PATH variable)
* bedtools (should be available via PATH variable)
* hisat2
* any of the supported trimming tools (trimmomatic, ...)
* any of the supported deduplication tools (fastuniq, ...)

<a name="venv"></a>
### Setting up virtual environment (recommended)
It is highly recommended to install package using in a virtual environment. This can be achieved via python's venv or conda envs.<br> Moreover, all of the required bioinformatic tools can be installed via conda as well (but mind the problems with [conda samtools installation](https://www.biostars.org/p/455593/)!) 


<a name="install1"></a>
### Clone and install package
clone repository from github and install:

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

### Config fields

* "rna_ids" : non-empty array of strings. Must contain names of files with RNA parts of reads, without directory names or file extensions. For example, for the file named <i>SRR9201799.fastq.gz</i> id is SRR9201799.
* "dna_ids" : non-empty array of strings.  Must contain names of files with DNA parts of reads, synchronized with rna_ids array.
* "base_dir" : string. Optional. Directory to run analysis in (one or several temporal directories will be created in this directory). Path should exist if provided. If not provided, current working directory will be taken.
* "input_dir" : string. Required. Path to directory with input DNA and RNA files.
* "output_dir" : string. Required. Path to directory to store results in. Will be created if doesnt exist.
* "cpus" : int. Number of threads to use.


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
        "params_list": [],
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

