import logging

import pandas as pd

GTF_NAMES = ('chr', 'type', 'start', 'end', 'strand', 'attrs')
logger = logging.getLogger(__name__)

class StrandCalc:
    def __init__(self,
                 input_dir: str,
                 gtf_annotation: str,
                 genes: str,
                 output_dir: str,
                 prefix: str) -> None:
        self.input_dir = input_dir
        self.output_dir = output_dir
        self.prefix = prefix
        self._load_genes(genes, gtf_annotation)
    
    def _load_genes(self,
                    genes: str,
                    gtf_annotation: str) -> None:
        with open(genes, 'r') as f:
            # will search for gene_name "DDX11L2"; OR gene_id "ENSG00000290825.1"; for all listed genes
            genes_list: str = '|'.join([f'"{line.strip()}";' for line in f])
        gene_annot = pd.read_csv(
            gtf_annotation, sep='\t', header=None, skiprows=5,
            usecols=[0,2,3,4,6,8], names=GTF_NAMES
        )
        gene_annot = gene_annot[gene_annot['type'] == 'gene']
        gene_annot = gene_annot[gene_annot['attrs'].str.contains(genes_list)]
        logger.debug(f'{gene_annot.shape[0]} genes selected from annotation file.')
        #ribgenes[8] = ribgenes[8].apply(lambda x: x.split('gene_name')[1].split(';')[0].split('"')[1])
        #ribgenes.set_index(8, inplace=True)
        self.gene_annot: pd.DataFrame = gene_annot
            

    def run(self) -> None:
        # https://stackoverflow.com/questions/17322109/get-dataframe-row-count-based-on-conditions
        pass
    '''
    for name in names:
        dat = pd.read_csv(f'~/data/imargi/bed/{name}.rna.bed', sep='\t', header=None)
        print(name, dat.shape[0])
        dat[0] = dat[0].apply(lambda x: hg38_dict[x])
        for item in ribgenes.index:
            same = dat[(dat[0]==ribgenes.at[item,0]) & (dat[5]==ribgenes.at[item,6]) & (dat[1] < ribgenes.at[item,4]) & (dat[2] > ribgenes.at[item,3])].shape[0]
            anti = dat[(dat[0]==ribgenes.at[item,0]) & (dat[5]!=ribgenes.at[item,6]) & (dat[1] < ribgenes.at[item,4]) & (dat[2] > ribgenes.at[item,3])].shape[0]
            res.at['imargi_'+name,item] = [same, anti]
    '''
    # save several otputs
    # save images (png and pdf)