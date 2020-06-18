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


prs.add_argument('-af','--add_filter',
                 required = False,
                 type = str,
                 nargs = 2,
                 default = [None,None],
                 help = ' '.join(['Additional filtering',
                                  'Column as first argument',
                                  'value as second',
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

if args.add_filter[0] is not None:
    add_filter = ds.ca[args.add_filter[0]] == args.add_filter[1]
else:
    add_filter = np.ones(ds.ca[label].shape[0]).astype(np.bool)


uni_types = np.unique(ds.ca[label].flatten())

use_cells = []

for type in uni_types:
    zidx = np.where((ds.ca[label].flatten() == type) & add_filter)[0]
    nz = zidx.shape[0]

    if nz < lower_bound:
        print(f"{type} | was discarded due to insufficient number of cells")
        continue

    elif nz >= upper_bound:
        zidx = np.random.choice(zidx,
                                size = upper_bound,
                                replace = False,
                                )

    fromz = zidx.shape[0]
    use_cells +=  zidx.tolist()

    print(f"{type} | Used {fromz} cells ")

use_cells = np.array(use_cells)
use_cells = np.sort(use_cells)


new_cnt = ds[:,use_cells].T
new_lbl = ds.ca[label].flatten()[use_cells]
new_barcodes = ds.ca['CellID'][use_cells]
_, idx = np.unique(new_barcodes,return_index = True)
new_cnt = new_cnt[idx,:]
new_lbl = new_lbl[idx]
new_barcodes = new_barcodes[idx]
genes = ds.ra['Gene']

new_cnt = pd.DataFrame(new_cnt,
                       index = pd.Index(new_barcodes),
                       columns = pd.Index(genes),
                      )

new_mta = pd.DataFrame(new_lbl,
                       index = pd.Index(new_barcodes),
                       columns = pd.Index(['bio_celltype']),
                       )


name,counts = np.unique(new_lbl,return_counts = True)

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
