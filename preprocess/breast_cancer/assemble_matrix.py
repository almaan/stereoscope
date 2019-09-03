#!/usr/bin/env python3

import pandas as pd
import numpy as np

import os.path as osp
from scipy.io import mmread

import argparse as arp

def stripper(x):
    if '.' in x:
        y = x.split('.')
        y = '.'.join(y[0:-1])
        return y
    else:
        return x

prs = arp.ArgumentParser()

prs.add_argument('-mp','--main_path',
                 type = str,
                 required = True,
                 help = '',
                )

prs.add_argument('-n','--n_cells',
                 type = int,
                 required = False,
                 default = 500,
                 help = '',
                )

prs.add_argument('-t','--tag',
                 type = str,
                 required = True,
                 help = '',
                )

prs.add_argument('-o','--output_dir',
                 type = str,
                 default = None,
                 required = False,
                 help = '',
                )


args = prs.parse_args()


tag = args.tag
pick_n_cells = args.n_cells

main_pth = args.main_path
output_dir = (args.main_path if args.output_dir is \
              None else args.output_dir)

mat_pth = osp.join(main_pth,"matrix.mtx")
brc_pth = osp.join(main_pth,"barcodes.tsv")
gen_pth = osp.join(main_pth,"genes.tsv")

genes = pd.read_csv(gen_pth, sep = '\t', index_col = 0, header = None,)
genes = genes.iloc[:,0].values.tolist()
barcodes = pd.read_csv(brc_pth, sep ='\t', index_col = None, header = None)
barcodes = barcodes.iloc[:,0].values.tolist()
barcodes = list(map(lambda x: tag + '_' + x,barcodes))

print(f" barcode length : {len(barcodes)}")
print(f" genes length : {len(genes)}")


mat = mmread(mat_pth).toarray().T
genesums = mat.sum(axis = 0)
cellsums = mat.sum(axis = 1)
cell_q = np.quantile(cellsums,0.1)
keep_cells = cellsums > cell_q
print(f" matrix size : {mat.shape}")
mat = mat[keep_cells,:]

genes = pd.Index(list(map(stripper,genes)))
uniq = (genes.duplicated() == False)
mat = mat[:,uniq]
genes = genes[uniq]
print(f" matrix size : {mat.shape}")

total_cells = mat.shape[0]
pick_n_cells = np.min((pick_n_cells,total_cells))
np.random.seed(1337)
picked_cells = np.random.choice(np.arange(total_cells), size = pick_n_cells, replace = False)
mat = mat[picked_cells,:]
barcodes = [ barcodes[x] for x in picked_cells ]




mat = pd.DataFrame(mat,
                   columns = genes,
                   index = barcodes,
                  )

mta = pd.DataFrame([tag] * mat.shape[0],
                   columns = ['bio_celltype'],
                   index = barcodes,
                  )

opth_cnt = osp.join(output_dir, '.'.join([tag,'cnt','processed','tsv']))
opth_mta = osp.join(output_dir, '.'.join([tag,'mta','processed','tsv']))

mat.to_csv(opth_cnt,
           sep = '\t',
           header = True,
           index = True,
          )

mta.to_csv(opth_mta,
           sep = '\t',
           header = True,
           index = True,
          )
