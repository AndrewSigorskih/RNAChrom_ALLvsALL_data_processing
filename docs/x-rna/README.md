<a name="head"></a>
# infer-xrna

This program tries to infer novel transcripts that do not correspond to any gene from user-provided annotation. It heavily relies on the results from 2 previous tools: *rnachromprocessing* and *detect-strand*.

## Table of Contents
1. [Procedure overview](#procedure)
2. [Inputs](#inputs)
3. [Usage](#usage)
4. [Config contents](#config)
5. [Outputs](#outputs)

<a name="inputs"></a>
## Inputs

* List of ids corresponding to RNA files to be processed.
* Reference gene annotation in **GTF** format. Reads that overlap with *any* interval from this annotation will be filtered out and not used for transcript inference.
* Main output of the [detect-strand](../detect-strand/README.md) tool from this package. This file should contain information about files to be processed.
* A path to a directory containing *bed* files, preferably obtained by using the [rnachromprocessing](../rnachromprocessing/README.md) tool from this package. Ids of the files should match with those listed as input ids.
* A path to a directory containing *fastq* files (plain text or gz-compressed). Ids of the files should match with those listed as input ids. Files should be preprocessed (preferably by using the [rnachromprocessing](../rnachromprocessing/README.md) tool from this package) and ready to be aligned to a reference genome.

<a name="usage"></a>
## Usage

```bash
infer-xrna [-h] -c CONFIG [-v/-vv]
```

Supported arguments:
```
-h, --help          show help message and exit
-c CONFIG, --config CONFIG
                    Configuration file in json or yaml format (required).
-v, --verbose       Verbosity level. By default little to none information is printed.
                    Use -v once to increase information logs about each step, and -vv to 
                    print every command that is being run.
```

<a name="config"></a>
## Config contents

Field name|Required|Description
---|---|---
rna_ids|Yes|non-empty array of strings. **Must** contain base names of files with RNA parts of reads, without directory names or file extensions. For example, for the file named <i>SRR9201799_1.fastq.gz</i> id is SRR9201799_1.
bed_input_dir|Yes|A path to a directory containing bed input files.
fq_input_dir|Yes|A path to a directory containing fastq input files.
output_dir|Yes|Path to directory to store results in. Will be created if doesn't exist.
base_dir|No|string. Directory to run analysis in (one or several temporal directories will be created in this directory). Path should exist if provided. If not provided, current working directory will be used instead.
cpus|No|int > 0. Number of tasks to run simultaneously on each step. Default is 1.
ouputs_prefix|No|string. Prefix for main output files and plots. Default value is "xrna".
keep_extras|No|non-empty array of strings. Defines a list of intermediate steps to be saved in output_dir after tool execution. Supported values:<br>- "ids": lists of ids of reads selected for further evaluation, one per input file.<br>- "fastq": input fastq files, filtered out by read ids from previous step and reverse-complemented if needed.<br>- "raw_bam": Bam files containing selected reads mapped to reference genome.<br>- "merged_bam": bam files from previous step, merged by biological replicas (defined by detect-strand output file) and sorted.<br>- "stringtie_raw": raw gtf files obtained by stringtie.<br>- "stringtie_merge": result of stringtie-merge procedure applied to files from previous step.<br>- "stringtie_cov": GTF files containing xrnas coverage (mean per-base, FPKM and TPM), one per merged bam file.
annotation|Yes|A sub-config containing paths to input annotation files (see below).
hisat|Yes|A sub-config containing information required to run hisat2 tool (see below).
stringtie|No|A sub-config containing information required to run stringtie tool (see below).

#### annotation sub-config:

Field name|Required|Description
---|---|---
gtf_annotation|Yes|Path to reference gene annotation in **GTF** format.
strand_info|Yes|Path to main output of the detect-strand tool.

#### hisat sub-config:

Field name|Required|Description
---|---|---
genome|Yes|Path ending with prefix of genome index files (will be passed to hisat's `-x` option; see more in [reference manual](http://daehwankimlab.github.io/hisat2/manual/#:~:text=number%3E%7D%20%5B%2DS%20%3Chit%3E%5D-,Main%20arguments,-%2Dx%20%3Chisat2%2Didx))
known_splice|Yes|Path to known-splice file (will be passed to hisat's `--known-splicesite-infile` option).
tool_path|No|Path to hisat2 executable if not available via PATH.
cpus|No|int > 0. Number of tasks to run simultaneously for alignment step. If not set, the "global" value will be used.
hisat_threads|No|int > 0. How many threads hisat will use **for every** simultaneously-run task. (Value to pass to hisat's `-p` option). Default is 1.

#### stringtie sub-config:
Field name|Required|Description
---|---|---
tool_path|No|Path to stringtie executable if not available via PATH.
cpus|No|int > 0. Number of tasks to run simultaneously for all parallel stringtie calls. If not set, the "global" value will be used.
stringtie_threads|No|int > 0. How many threads stringtie will use **for every** simultaneously-run task. (Value to pass to stringtie's `-p` option). Default is 1.

[Back to top](#head)
