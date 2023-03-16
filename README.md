# RNAChrom_ALLvsALL_data_processing

Pachage for processing data of several "ALL-vs-ALL" RNA-Chromatin interactions capturing experiments, e.g. RADICL, GRID, iMARGI, CHAR, etc.

## Table of Contents
1. [Installation](#installation)
    1. [Required software](#software1)
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

### Setting up virtual environment (recommended)
It is highly recommended to install package using in a virtual environment. This can be achieved via python's venv or conda envs.<br> Moreover, all of the required bioinformatic tools can be installed via conda as well (but mind the problems with [conda samtools installation](https://www.biostars.org/p/455593/)!) 



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

