#!/usr/bin/env python3

import os.path as osp
import os
import datetime
import sys
import re

import numpy as np
import pandas as pd
import loompy as lp

import argparse as arp

prs = arp.ArgumentParser()

prs.add_argument('-lf','--loompy_file',
                 required = True,
                 type = str,
                 help = ' '.join(['Full path to loompy',
                                  'file.',
                                 ]))

prs.add_argument('-o','--output_dir',
                 required = False,
                 type = str,
                 default = '.',
                 help = ' '.join(['Path to output directory.',
                                  'If not specified curret work',
                                  'directory is used',
                                 ]))

prs.add_argument('-lb','--lower_bound',
                 required = False,
                 type = int,
                 default = 0,
                 help = ' '.join(['Lower bound for number',
                                  'of cells to be present',
                                  'in a cell type.',
                                  'Types with less than specified',
                                  'value will be omitted from set',
                                 ]))

prs.add_argument('-ub','--upper_bound',
                 required = False,
                 type = int,
                 default = None,
                 help = ' '.join(['Upper bound for number',
                                  'of cells to be included',
                                  'from a cell type. Types'
                                  'with more cells than specified'
                                  'value will be randomly',
                                  'subsampled',
                                 ]))

prs.add_argument('-cn','--column_name',
                 required = False,
                 type = str,
                 default = 'Subtype',
                 help = ' '.join(['Name of column in loompy',
                                  'file which specifies',
                                  'type label',
                                 ]))

args = prs.parse_args()

np.random.seed(1337)

loom_pth = args.loompy_file
lower_bound = args.lower_bound
upper_bound = (args.upper_bound if args.upper_bound is not None else np.inf)
out_dir = args.output_dir

if not osp.exists(out_dir):
    os.mkdir(out_dir)

label = args.column_name

tag = str(datetime.datetime.now())
tag = re.sub('-| |:|\\.','',tag)

print(f"Unique Identifier for set >> {tag}")

try:
    ds = lp.connect(loom_pth)
except OSError:
    print(f'The file {loom_pth} either does not exist or',
          f" you might already have established a connection",
          f" to the file.")
    sys.exit(-1)

genes = ds.ra['Gene']
ngenes = genes.shape[0]

uni_types = np.unique(ds.ca[label])
n_uni_types = uni_types.shape[0]

cell_list = []
barcode_list = []
idx_list = [0]
use_labels = []

for z in range(n_uni_types):
    type = uni_types[z]
    zidx = np.where(ds.ca[label] == type)[0]
    nz = zidx.shape[0]

    if nz < lower_bound:
        print(f"{type} | was discarded due to insufficient number of cells")
        continue

    elif nz > upper_bound:
        zidx = np.random.choice(np.where(zidx)[0],
                                    size = upper_bound,
                                    replace = False,
                                   )

    zidx.sort()
    fromz = zidx.shape[0]

    use_labels += [type] * fromz
    cell_list.append(ds[:,zidx])
    idx_list.append(fromz)

    barcode_list += ds.ca['CellID'][zidx].tolist()

    print(f"{type} | Used {fromz} cells ")


idx_list = np.array(idx_list)
pos  = np.cumsum(idx_list)
n_uni_types = len(cell_list)

new_cnt = np.zeros((pos[-1],ngenes))
print(new_cnt.shape)

for z in range(n_uni_types - 1):
    new_cnt[pos[z]:pos[z + 1],:] = cell_list[z].T

new_cnt = pd.DataFrame(new_cnt,
                       index = pd.Index(barcode_list),
                       columns = pd.Index(genes),
                      )

new_mta = pd.DataFrame(use_labels,
                       index = pd.Index(barcode_list),
                       columns = pd.Index(['type']),
                       )


name,counts = np.unique(use_labels,return_counts = True)

stats = pd.DataFrame(counts,
                     index = pd.Index(name),
                     columns = pd.Index(['members']),
                     )

print(stats)
print(f"Assembled Count Matrix with dimension >> nrow : {new_cnt.shape[0]} | ncol :  {new_cnt.shape[1]}")
print(f"Assembled Meta File with >> nrow : {new_mta.shape[0]}")

opth_cnt = osp.join(out_dir,
                    '.'.join([tag,
                             'cnt_data.tsv',
                            ]
                           )
                   )


opth_mta = osp.join(out_dir,
                    '.'.join([tag,
                             'mta_data.tsv',
                            ]
                           )
                   )

opth_sts = osp.join(out_dir,
                    '.'.join([tag,
                             'stats.tsv',
                            ]
                           )
                   )


new_cnt.to_csv(opth_cnt,
               sep = '\t',
               header = True,
               index = True,
               index_label = 'cell',
              )

new_mta.to_csv(opth_mta,
               sep = '\t',
               header = True,
               index = True,
               index_label = 'cell',
              )

stats.to_csv(opth_sts,
             sep = '\t',
             header = True,
             index = True,
             index_label = 'cell',
              )
ds.close()
