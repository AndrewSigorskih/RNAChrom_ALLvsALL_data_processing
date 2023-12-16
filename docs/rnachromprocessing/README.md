<a name="head"></a>
# rnachromprocessing

This program performs several consecutive steps of data processing as follows:

1. Restriction sites handling. Filtering reads by starting bases and/or finishing the fragmented restriction site. Several "procedure presets" are available for this step based on previous papers, however users can also use their own scripts by passing those as innputs (see "Config contents" section for more information).
2. Dedupication step. PCR duplicates are identified and removed.
3. Trimming. Trimming of the input fastq files by quality. 
4. Mapping RNA and DNA reads to reference genome separately using the chosen alignment software
5. Filtering resulting bam files: only reads that were mapped with N (user-defined; default 2) or fewer mismatches will pass this stage.
6. Filtering resulting bam files: only reads that were uniquelly mapped will pass this stage. Filtered bam files are converted to bed format.
7. Joining surviving read pairs into resulting contacts table.

Each step from 1 to 3 can be skipped. The pipeine interface provides opportunity to fine-tune the behaviour for almost each step, hovewer in case of running with defaults the user's input can be minimalistic.


## Table of Contents
1. [Usage](#processingusage)
2. [Config contents](#processingconfig)
    * [Providing custom scripts for pipeline stages](#custom)
3. [Minimal config example](#minimalconfig)
4. [Advanced config example](#advancedconfig)
5. [Running only specified stage](#singlestage)
6. [Outputs](#processingoutputs)

<a name="processingusage"></a>
### Usage

```
rnachromprocessing [-h] -c CONFIG [-v/-vv] [-s STAGE] [--input_dir INPUT_DIR] [--output_dir OUTPUT_DIR]
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
                    Use -v once to get information logs about each step, and -vv to 
                    print detailed information about every command that is being run.
--input_dir INPUT_DIR
                    Specify input directory. Overrides the "input_dir" field in config.
--output_dir OUTPUT_DIR
                    Specify output directory. Overrides the "output_dir" field in config.    
```

All inputs and parameters are passed using the </i>config</i> argument -- a path to existing file in <b>yaml</b> or <b>json</b> format.

Config should contain several required fields and may contain several additional fields for fine-tuning. Complete config specification and examples are listed below.

**Config examples can be found in `configs/rnachromprocessing` directory of this repository.**

<a name="processingconfig"></a>
### Config contents

Input config is expected to be in a form of yaml/json file with several "global" fields and several sub-configs corresponding to individual stages of the pipeline. Detailed explanation of each parameter is listed in the tables below:

#### "Global" options:

Field name|Required|Description
---|---|---
rna_ids|Yes|non-empty array of strings. **Must** contain base names of files with RNA parts of reads, without directory names or file extensions. For example, for the file named <i>SRR9201799_1.fastq.gz</i> id is SRR9201799_1.
dna_ids|Yes|non-empty array of strings.  **Must** contain names of files with DNA parts of reads, synchronized with rna_ids array.
base_dir|No|string. Directory to run analysis in (one or several temporal directories will be created in this directory). Path should exist if provided. If not provided, current working directory will be used instead.
input_dir|Yes|Path to directory with input DNA and RNA files. Should exits.
output_dir|Yes|Path to directory to store results in. Will be created if doesn't exist.
cpus|No|int > 0. Number of tasks to run simultaneously on each step. Default is 1.
keep|No|non-empty array of strings. Should contain names of steps that will have their results saved into "output_dir". Supported values (can be listed in any order): ["rsites", "dedup", "trim", "align", "bam", "bed", "contacts"]. By default only "trim", "bed" and "contacts" steps are saved.
rsites|No|A sub-config containing parameters for the "rsites" stage (see below).
dedup|No|A sub-config containing parameters for the "dedup" stage (see below).
trim|No|A sub-config containing parameters for the "trim" stage (see below).
align|Yes|A sub-config containing parameters for the "align" stage (see below).
bam|No|A sub-config containing parameters for the "bam" stage (see below).
bed|No|A sub-config containing parameters for the "bed" stage (see below).
contacts|No|A sub-config containing parameters for the "contacts" stage (see below).
stats|No|A sub-config containing parameters for statistics calculation (see below).


#### rsites stage subconfig:

A sub-config for restriction-site handling step. Supported options:


Field name|Required|Description
---|---|---
cpus|No|int > 0. Number of tasks to run simultaneously for this step. If not set, the "global" value will be used.
type|No|What action to perform on read pairs during this step. Default is "skip". Supported values:<br>- "imargi" : Save read pair if DNA read starts with CT or NT and remove first 2 bases from RNA reads. Reads in files should be synchronized.<br>- "grid" : Save read pair if DNA read ends with **rsite_bgn** bases (default AG). Add **rsite_end** bases (default CT) to the end of selected DNA reads (quality is copied from **rsite_bgn** bases). Reads in files should be synchronized.<br>- "custom" : Use custom script provided by user. See [custom script section](#custom) for more information.<br>- "skip" : do not perform any data transformations during this step.
rsite_bgn|No|Beginning of restriction site that will be used in "grid"-like processing. Default "AG".
rsite_end|No|Missing end of restriction site that will be used in "grid"-like processing. Default "CT".
tool_path|No|Will be only considered if "custom" type of action is chosen for "type" option. Path to user's script that will manage all operations on restriction sites for read pairs.

#### dedup stage subconfig:

A sub-config for read deduplication step. Supported options:

Field name|Required|Description
---|---|---
cpus|No|int > 0. Number of tasks to run simultaneously for this step. If not set, the "global" value will be used.
tool|No|What deduplication tool to use. Default is "fastuniq". Supported options:<br>- "fastuniq" : use fastuniq tool.<br>- "fastq-dupaway": use our memory-efficient deduplication tool "fastq-dupaway".<br>- "custom" : Use custom script provided by user. See [custom script section](#custom) for more information.<br>- "skip" : do not perform deduplication step.
tool_path|No|Path to tool executable. Use if your tool is not available via PATH. If this option is not set, package will try to search for executable location using the `shutil.which()` function with the name of tool provided in previous option. If "custom" tool is chosen, path to your script must be provided here.
tool_params|No|Sub-config with parameters to be passed to deduplication tool(s). Supported fields:<br>- "memlimit": int in range [500, 10240] (default 4096). Supported by fastq-dupaway tool. Memory limit in megabytes for every call of fastq-dupaway to use.<br>- "comparison": either "tight" or "loose" (default). Choose sequance comparison logic for fastq-dupaway tool.

#### trim stage subconfig:

A sub-config for read trimming step, should be found under the name "trim". Supported options:

Field name|Required|Description
---|---|---
cpus|No|int > 0. Number of tasks to run simultaneously for this step. If not set, the "global" value will be used.
tool|No|What trimming tool to use. Default is "trimmomatic". Supported options:<br>- "trimmomatic"<br>- "custom" : Use custom script provided by user. See [custom script section](#custom) for more information.<br>- "skip" : do not perform trimming step.
tool_path|No|Path to tool executable. Use if your tool is not available via PATH. If this option is not set, package will try to search for executable location using the `shutil.which()` function with the name of tool provided in previous option. If "custom" tool is chosen, path to your script must be provided here.
tool_params|No|Tool-specific parameters in form of a sub-config. For "Trimmomatic" tool supported parameters are "window" (default 5), "qual_th" (default 26) and "minlen" (default 15).

#### align stage subconfig:

A sub-config for read alignment step, should be found under the name "align". Supported options:

Field name|Required|Description
---|---|---
cpus|No|int > 0. Number of tasks to run simultaneously for alignment step. If not set, the "global" value will be used.
tool|No|What alignment tool to use. Default is "hisat2". Supported options:<br>- "hisat2": use hisat2 aligner with `--no-spliced-alignment` mode for DNA reads and `--dta-cufflinks` for RNA reads.<br>- "bwa": [not implemented yet]<br>- "star": [not implemented yet]<br>- "custom": Use custom script provided by user. See [custom script section](#custom) for more information.
tool_path|No|Path to tool executable. Use if your tool is not available via PATH. If this option is not set, package will try to search for executable location using the `shutil.which()` function with the name of tool provided in previous option. If "custom" tool is chosen, path to your script must be provided here.
dna_genome_path|**Yes**|Path ending with prefix of genome index files for DNA reads alignment.
rna_genome_path|No|Path ending with prefix of genome index files for RNA reads alignment. If not provided, the dna_genome_path will be used here as well.
known_splice|**Yes** for hisat2|Path to known-splice file (will be passed to hisat's `--known-splicesite-infile` option).
tool_threads|No|int > 0. How many threads aligner tool will use **for every** simultaneously-run task. Default is 1. Total maximum core usage for the "align" stage will be cpus * tool_threads.

#### bam stage subconfig:

A sub-config for bam-filtering step, should be found under the name "bam". Supported options:

Field name|Required|Description
---|---|---
cpus|No|int > 0. Number of tasks to run simultaneously for bam-filtering step. If not set, the "global" value will be used.
max_mismatch|No|int > 0. Maximum number of mismatches allowed in bam record to pass this step. Default value is 2.

#### bed stage subconfig:

A sub-config for bam-to-bed conversion step, should be found under the name "bed". Supported options:

Field name|Required|Description
---|---|---
cpus|No|int > 0. Number of tasks to run simultaneously for bam-to-bed conversion step. If not set, the "global" value will be used.

#### contacts stage subconfig:

A sub-config for contacts building step, should be found under the name "contacts". Supported options:

Field name|Required|Description
---|---|---
cpus|No|int > 0. Number of tasks to run simultaneously for contacts-building step. If not set, the "global" value will be used.
mode|No|memory-efficiency toggling option. Default is "fast". Supported values:<br>- "fast" : load all data in ram at once.<br>- "low-mem" : process data in a slower fashion, but more RAM-friendly.


#### Statistic accumulation stage subconfig:

A sub-config for controling surviving reads statistics calculation, should be found under the name "stats". Supported options:

Field name|Required|Description
---|---|---
cpus|No|int > 0. Number of tasks to run simultaneously for every stats-calculating step. If not set, the "global" value will be used.
mode|No|Whether or not count surviving reads statistics. Default is "default". Supported values:<br>- "default" : gather sirviving reads counts after each step, write resulting table in the output folder.<br>- "skip" : do not perform this stage.
prefix|No|Name of the resulting table (without the .tsv prefix) that will be written in the output folder. Default is "stats".


<a name="custom"></a>
### Providing custom scripts for pipeline stages.

Several pipeline stages (rsites, dedup, trim, align) may accept custom user script instead of running pre-defined behaviour. Such script should comply with the following rules:

* It should be an executable that accepts exactly 4 positional arguments: input_dna_file, input_rna_file, output_dna_file, output_rna_file.
During pipeline execution, it will be run as follows:

```<script_name> input_dna_file input_rna_file output_dna_file output_rna_file```

* Input and output files for fastq-related stages should have the same compression status, i.e. if the input data was gzip-compressed, output should be compressed as well, and vice versa.

* Output files for the "align" stage should be in bam format **without unaligned reads**.

* As the pipeline will most likely execute your script in parallel, make sure that temporary files (if any) do not overlap and are located in the current stage temporary directory: 
    * Don't use generic names as "file.tmp", instead use the input/output file paths to create names for temporary files, e.g. `${OUTPUT_DNA}.tmp`
    * Make sure to remove all such temporary files at the end of script execution.


<a name="minimalconfig"></a>
### Minimal config example

<details>

<summary>An example of config file with minimalistic settings in yaml format:</summary>

```
rna_ids:
- SRR9201799
- SRR9201801
- SRR9201803
- SRR9201805
- SRR9201807
dna_ids:
- SRR9201800
- SRR9201802
- SRR9201804
- SRR9201806
- SRR9201808
input_dir: /gpfs/asigorskikh/data/radicl/gz
output_dir: new_results
cpus: 16
align:
    cpus: 2
    tool: hisat2
    dna_genome_path: /gpfs/asigorskikh/data/genomes/mm10/mm10
    known_splice: /gpfs/asigorskikh/data/genes/gencode.vM32.ss
    tool_threads: 8
stats:
    cpus: 8
    prefix: radicl
    mode: default

```

</details>

<a name="advancedconfig"></a>
### Advanced config example

<details>

<summary>An example of heavily  customized config file in yaml format:</summary>

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
        "cpus": 2
    },
    "rsites": {
        "type": "imargi"
    },
    "trim": {
        "tool" : "trimmomatic",
        "tool_path" : "path/to/trimmomatic/jar",
        "tool_params": {
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

</details>

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

* **output_dir/stats::prefix**.tsv: Table with passing reads statistic.


[Back to top](#head)
