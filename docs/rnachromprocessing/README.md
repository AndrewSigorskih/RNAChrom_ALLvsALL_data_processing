<a name="head"></a>
# rnachromprocessing

This program performs several consecutive steps of data processing as follows:

1. Restriction sites handling. Filtering reads by starting bases and/or finishing the fragmented restriction site. Several "procedure presets" are available for this step based on previous papers, however users can also use their own scripts by passing those as innputs (see "Config contents" section for more information).
2. Dedupication step. PCR duplicates are identified and removed.
3. Trimming. Trimming of the input fastq files by quality. The default behaviour for this step is using Trimmomatic with parameters **SLIDINGWINDOW:5:26 MINLEN:15**.
4. Mapping RNA and DNA reads to reference genome separately using the Hisat2 software:
    * **--no-spliced-alignment -k 100 --no-softclip** parameters for DNA read parts,
    * **--dta-cufflinks -k 100 --no-softclip --known-splicesite-infile** for RNA read parts.
5. Filtering resulting bam files: only reads that were uniquely mapped with 2 or fewer mismatches will pass this stage.
6. Converting filtered bam files to bed and joining them into a resulting contacts table.
7. Overall statistics calculation: how many pairs of reads survived each step.

Each step from 1 to 3 can be skipped. The pipeine interface provides opportunity to fine-tune the behaviour for almost each step, hovewer in case of running with defaults the user's input can be minimalistic.


## Table of Contents
1. [Usage](#processingusage)
2. [Config contents](#processingconfig)
3. [Minimal config example](#minimalconfig)
4. [Advanced config example](#advancedconfig)
5. [Running only specified stage](#singlestage)
6. [Outputs](#processingoutputs)

<a name="processingusage"></a>
### Usage

```
rnachromprocessing [-h] -c CONFIG [-v/-vv] [-s STAGE]
```

Supported arguments:
```
-h, --help          show help message and exit
-c CONFIG, --config CONFIG
                    Configuration file in json format (required).
-s STAGE, --stage STAGE
                    Run only specified stage of the pipeline.
                    Supported stages: {dedup, rsites, trim, hisat, bam, bed, contacts}.
-v, --verbose       Verbosity level. By default little to none information is printed.
                    Use -v once to increase information logs about each step, and -vv to 
                    print every command that is being run.
```

All inputs and parameters are passed using the </i>config</i> argument -- a path to existing file in a <b>.json</b> format.

Config should contain several required fields and may contain several additional fields for fine-tuning. Complete config specification and examples are listed below.

**Config examples can be found in `configs/rnachromprocessing` directory.**

<a name="processingconfig"></a>
### Config contents

Input config is expected to be in a form of json file with several "global" fields and several sub-configs corresponding to individual stages of the pipeline. Detailed explanation of each parameter is listed in the tables below:

#### "Global" options:

Field name|Required|Description
---|---|---
rna_ids|Yes|non-empty array of strings. **Must** contain names of files with RNA parts of reads, without directory names or file extensions. For example, for the file named <i>SRR9201799_1.fastq.gz</i> id is SRR9201799_1.
dna_ids|Yes|non-empty array of strings.  **Must** contain names of files with DNA parts of reads, synchronized with rna_ids array.
base_dir|No|string. Directory to run analysis in (one or several temporal directories will be created in this directory). Path should exist if provided. If not provided, current working directory will be used instead.
input_dir|Yes|Path to directory with input DNA and RNA files.
output_dir|Yes|Path to directory to store results in. Will be created if doesn't exist.
cpus|No|int > 0. Number of tasks to run simultaneously on each step. Default is 1.
keep|No|non-empty array of strings. Should contain names of steps that will have their results saved into "output_dir". Supported values (can be listed in any order): ["rsites", "dedup", "trim", "hisat", "bam", "bed", "contacts"]. By default only "bam" and "contacts" steps are saved.
stats|No|either "skip", "default" or "full". Optional. If "skip" value is chosen no passing reads statistic will be collected. Default behaviour is calculating number of **pairs** of reads that passed rsite, dedup and trim steps as well as resulting number of contacts. "full" mode will also calculate number of pairs of reads that were mapped  to reference genome regardless of any filtration and number of mapped reads that had 2 or less mismatches, however running this mode will take a lot of extra time.
stats_prefix|No|string. Optional, default value is "stats". Prefix name for output table with passing reads statistic.


#### rsites stage subconfig:

A sub-config for restriction-site handling step, should be found under the name "rsites". Supported options:


Field name|Required|Description
---|---|---
cpus|No|int > 0. Number of tasks to run simultaneously for this step. If not set, the "global" value will be used.
type|No|What action to perform on read pairs during this step. Default is "skip". Supported values:<br>- "imargi" : Save read pair if DNA read starts with CT or NT and remove first 2 bases from RNA reads. Reads in files should be synchronized.<br>- "grid" : Save read pair if DNA read ends with **rsite_bgn** bases (default AG). Add **rsite_end** bases (default CT) to the end of selected DNA reads (quality is copied from **rsite_bgn** bases). Reads in files should be synchronized.<br>- "custom" : Use custom script provided by user.<br>- "skip" : do not perform this step.
rsite_bgn|No|Beginning of restriction site that will be used in "grid"-like processing. Default "AG".
rsite_end|No|Missing end of restriction site that will be used in "grid"-like processing. Default "CT".
tool_path|No|Will be only considered if "custom" type of action is chosen for "type" option. Path to user's script that will manage all operations on restriction sites for read pairs. Script should be executable and take exactly 4 positional arguments: `dna-input-file`, `rna-input-file`, `dna-output-file`, `rna-output-file`.

#### dedup stage subconfig:

A sub-config for read deduplication step, should be found under the name "dedup". Supported options:

Field name|Required|Description
---|---|---
tool|No|What deduplication tool to use. Default is "fastuniq". Supported options:<br>- "fastuniq" : use fastuniq tool.<br>- "fastq-dupaway": use our memory-efficient deduplication tool "fastq-dupaway".<br>- "skip" : do not perform deduplication step.
tool_path|No|Path to tool executable. Use if your tool is not available via PATH. If this option is not set, package will try to search for executable location using the `shutil.which()` function with the name of tool provided in previous option.
cpus|No| int > 0. Number of tasks to run simultaneously for deduplication step. If not set, the "global" value will be used.
params|No|Sub-config with parameters to be passed to deduplication tool(s). Supported fields:<br>- "memlimit": int in range [500, 10240]. Supported by fastq-dupaway tool. Memory limit in megabytes for every call of fastq-dupaway to use.

#### trim stage subconfig:

A sub-config for read trimming step, should be found under the name "trim". Supported options:

Field name|Required|Description
---|---|---
tool|No|What trimming tool to use. Default is "trimmomatic". Supported options:<br>- "trimmomatic"<br>- "skip" : do not perform trimming step.
tool_path|No|Path to tool executable. Use if your tool is not available via PATH. If this option is not set, package will try to search for executable location using the `shutil.which()` function with the name of tool provided in previous option.
cpus|No|int > 0. Number of tasks to run simultaneously for deduplication step. If not set, the "global" value will be used.
params|No|Tool-specific parameters in form of a sub-config. For "Trimmomatic" tool supported parameters are "window", "qual_th" and "minlen" (int > 0 each).

#### hisat stage subconfig:

A sub-config for read alignment step, should be found under the name "hisat". Supported options:

Field name|Required|Description
---|---|---
genome|Yes|Path ending with prefix of genome index files (will be passed to hisat's `-x` option; see more in [reference manual](http://daehwankimlab.github.io/hisat2/manual/#:~:text=number%3E%7D%20%5B%2DS%20%3Chit%3E%5D-,Main%20arguments,-%2Dx%20%3Chisat2%2Didx))
known_splice|Yes|Path to known-splice file (will be passed to hisat's `--known-splicesite-infile` option).
cpus|No|int > 0. Number of tasks to run simultaneously for deduplication step. If not set, the "global" value will be used.
hisat_threads|No|int > 0. How many threads hisat will use **for every** simultaneously-run task. (Value to pass to hisat's `-p` option). Default 1.

#### contacts stage subconfig:

A sub-config for contacts building step, should be found under the name "contacts". Supported options:

Field name|Required|Description
---|---|---
cpus|No|int > 0. Number of tasks to run simultaneously for deduplication step. If not set, the "global" value will be used.
mode|No|memory-efficiency toggling option. Default is "fast". Supported values:<br>- "fast" : load all data in ram at once.<br>- "low-mem" : process data in a slower fashion, but more RAM-friendly.


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
    "rna_ids": ["SRR12462453_1", "SRR8206679_1", "SRR8206680_1", "SRR9900121_1", "SRR9900122_1"],
    "dna_ids": ["SRR12462453_2", "SRR8206679_2", "SRR8206680_2", "SRR9900121_2", "SRR9900122_2"],
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
    },
    "contacts" : {
        "mode" : "low-mem"
    }
}
```

<a name="singlestage"></a>
### Running only specified stage of the pipeline

<i>rnachromporcessing</i> supports the ability to run each stage of the pipeline individually, which is achieved by running with `-s/--stage` option. The config file can be the same as the one used for "full" mode, and in that case all sub-configs for other stages will be simply ignored.

The `input_dir` field should reference folder with valid inputs for chosen stage.


<a name="processingoutputs"></a>
### Outputs

* **output_dir**/contacts/*.tab:
    Tables in tsv format named after RNA files. Columns:
    * 'id': read pair unique ID.
    * 'rna_chr', 'rna_bgn', 'rna_end', 'rna_strand', 'rna_cigar': RNA read information from corresponfing bed file.
    * 'dna_chr', 'dna_bgn', 'dna_end', 'dna_strand', 'dna_cigar': DNA read information from corresponfing bed file.

* **output_dir/stats_prefix**.tsv: Table with passing reads statistic.


[Back to top](#head)
