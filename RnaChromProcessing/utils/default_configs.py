dedup_default_cfg = {
    'tool': 'fastuniq',
    'params': [],
    'cpus': None
}

rsites_default_cfg = {
    'cpus': None,
    'type': 'skip'
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
    'known_splice': None
}