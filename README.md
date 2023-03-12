# RNAChrom_ALLvsALL_data_processing

Pachage for processing data of several "ALL-vs-ALL" RNA-Chromatin interactions capturing experiments, e.g. RADICL, GRID, iMARGI, CHAR, etc.

[!Work in progress]
## Installation

## Usage


### Minimal config example

```

```

### Advanced config example

```
{
    'rna_ids': [],
    'dna_ids': [],
    'schema': 'radicl',
    'base_dir' : '.',
    'input_dir': './fastq',
    'output_dir': 'results',
    'cpus': 5,
    'keep' : ['contacts', 'bam'],
    'dedup': {
        'tool': 'fastuniq',
        'tool_path': '/path/to/fastuniq',
        'params_list': [],
        'cpus': None
    },
    'rsites': {
        'cpus': None,
        'type': 'skip'
    },
    'trim': {
        'tool' : 'trimmomatic',
        'tool_path' : 'path/to/trimmomatic/jar',
        'cpus': None,
        'params': {
            'window': 5,
            'qual_th': 26,
            'minlen': 15
        }
    },
    'hisat': {
        'tool_path' : 'path/to/hisat',
        'genome': 'PATH/TO/GENOMEDIR',
        'known_splice': 'PATH/TO/KNOWN/SPLICE',
        'cpus': None
    }
}
```
