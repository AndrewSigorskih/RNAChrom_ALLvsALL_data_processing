<a name="head"></a>
# detect-strand

This program takes contacts files from produced py rnachromprocessing as input and checks whether orientation of RNA parts of contacts was inverted or lost during sequencing.

## Table of Contents
1. [Usage](#strandusage)
2. [Config contents](#strandconfig)
3. [Outputs](#outputs)

The following procedure is performed:
* User provides gene annotation file in GTF format and file with a list of selected genes. This list can contain gene names or gene ids (from 9th field of gene annotation file).
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
-v, --verbose       Verbosity level. By default little to none information is printed.
                    Use -v once to increase information logs about each step, and -vv to 
                    print every command that is being run.
```

**Examples of configs and genes lists for human and mouse datasets can be found in `configs/strand` directory.**

<a name="strandconfig"></a>
### Config contents

Field name|Required|Description
---|---|---
input_dir|Yes|path to input directory with contacts tables.
output_dir|Yes|path to directory to store outputs in. Will be created if not exist.
gtf_annotation|Yes|path to gene annotation file in GTF format.
genes_list|Yes|path to a file with selected genes names (or ids), one per line.
prefix|No|prefix name for all output files. Default  is "strand".
cpus|No|int > 0. Number of tasks to run simultaneously. Default is 1.
plots_format|No|Either "png" (default) or "svg". Specifies format for saving resulting plots.
exp_groups|Yes|sub-config with lists of input file ids divided by groups (experiment types, for example). See config examples for more information.


<a name="outputs"></a>
### Outputs

* **output_dir**/**prefix**_wins.tsv:
    Main output. A table in .tsv format, containing information about wins and loses of "same" and "antisense" strands for each dataset.<br>
    This file is a required input for infer-xrna program.

* **output_dir**/**prefix**_raw_counts.tsv:
    Table in .tsv format, containing pairs of "same" and "anti" read counts for each gene for each dataset.<br>
    Columns represent genes, multiindex represents experimental groups derived from the "exp_groups" config field, each value is a pair of integers in a form of a pythonic tuple.<br>
    In order to read the table properly use the following code snippet:<br>
    `pd.read_csv('PATH/TO/PREFIX_raw_counts.tsv', index_col=(0, 1), sep='\t').applymap(eval)`

* **output_dir**/**prefix**_wins.png:
    A barplot that visualizes the "wins.tsv" table. Can be in either png or svg format. For each dataset, two bars are drawn: the bar with positive height represents the number of times the "same" strand won out of all genes, and the bar with negative height represents the number of times the "antisense" strand won. Bar colors represent experimental groups.

* **output_dir**/**prefix**_cov_boxplot.png:
    Boxplot that visualizes the "raw_counts.tsv" table. Can be in either png or svg format. For each dataset, two boxplots are drawn, representing the corresponding coverage of each gene by RNA parts on same and on antisense strand. Doverage values are log10-transformed.

[Back to top](#head)
