#!/usr/bin/env python3

import pandas as pd
import numpy as np

import os.path as osp
import re
import datetime
from scipy.io import mmread

import argparse as arp

tag = str(datetime.datetime.now())
tag = re.sub('-| |:|\\.','',tag)

prs = arp.ArgumentParser()

prs.add_argument('-cp','--count_paths',
                 type = str,
                 nargs = '+',
                 required = True,
                 help = '',
                )

prs.add_argument('-mp','--meta_paths',
                 type = str,
                 nargs = '+',
                 required = True,
                 help = '',
                )
prs.add_argument('-o','--output_dir',
                 type = str,
                 default = None,
                 help = '',
                )

args = prs.parse_args()

output_dir = (osp.dirname(args.count_paths[0]) if args.output_dir is \
              None else args.output_dir)

genes = pd.Index([])
barcodes = pd.Index([])
cnt_list = []

start_pos = [0]
new_mta = pd.DataFrame([])

for ii in range(len(args.count_paths)):

    ct = pd.read_csv(args.count_paths[ii],
                     sep = '\t',
                     header = 0,
                     index_col = 0,
                    )

    cnt_list.append(ct)

    mt = pd.read_csv(args.meta_paths[ii],
                     sep = '\t',
                     header = 0,
                     index_col = 0,
                    )

    new_mta = pd.concat([new_mta,mt])

    genes = genes.union(ct.columns)
    barcodes = barcodes.union(ct.index)
    start_pos.append(ct.shape[0])

start_pos = np.cumsum(np.array(start_pos))
new_cnt  = pd.DataFrame(np.zeros((start_pos[-1],genes.shape[0])),
                        columns = genes,
                        index = barcodes,
                       )

for t in range(len(cnt_list) - 1):
    start = start_pos[t]
    end = start_pos[t+1]
    inter = cnt_list[t].columns.intersection(genes)
    new_cnt.loc[start:end,inter] = cnt_list[t].loc[:,inter].values


opth_cnt = osp.join(output_dir,'.'.join([tag,'cnt','tsv']))
opth_mta = osp.join(output_dir,'.'.join([tag,'mta','tsv']))

new_cnt.to_csv(opth_cnt,
           sep = '\t',
           header = True,
           index = True,
          )

new_mta.to_csv(opth_mta,
           sep = '\t',
           header = True,
           index = True,
           )

