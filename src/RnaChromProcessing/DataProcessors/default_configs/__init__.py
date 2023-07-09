SUBDIR_LIST = ('dedup', 'rsites', 'trim', 'hisat', 'bam', 'bed', 'contacts')

dedup_default_cfg = {
    'tool': 'fastuniq',
    'params': [],
    'cpus': None
}

rsites_default_cfg = {
    'cpus': None,
    'type': 'skip',
    'rsite_bgn': 'AG',
    'rsite_end': 'CT'
}

trim_default_cfg = {
    'tool': 'trimmomatic',
    'params': {
        'window': 5,
        'qual_th': 26,
        'minlen': 15
    },
    'cpus': None
}

hisat_default_cfg = {
    'cpus': None,
    'genome': None,
    'known_splice': None,
    'hisat_threads': 1
}

contacts_default_cfg = {
    'cpus': None,
    'mode': 'fast',
    'chunksize': 1_000_000
}