#!/usr/bin/env python3

import os.path as osp
import os

import pandas as pd
import numpy as np

import argparse as arp

prs = arp.ArgumentParser()
prs.add_argument('-i','--input')
prs.add_argument('-at','--add_total',
                 type = str,
                 default = None,
                )
args = prs.parse_args()




tbl = pd.read_csv(args.input,
                  sep = '\t',
                  index_col = None,
                  header = 0,
                 )

if args.add_total is not None:
    tot =  np.sum(tbl.loc[:,args.add_total].values)
    tot_row = tbl.iloc[0,:]
    pos = tbl.columns == args.add_total
    tot_row.iloc[0] = 'Total'
    tot_row.iloc[pos] = str(tot)
    tbl = tbl.append(tot_row)
    print(tot_row)


print(f"Shape of Table {tbl.shape}")
n_cols = tbl.shape[1]
prms = {'column_format':'|l|' + 'c|'*(n_cols-1),
        'bold_rows':True,
        'index':False,
        }

tbl_ltx  = tbl.to_latex(buf = None,
                       **prms)

tbl_ltx = tbl_ltx.splitlines()
tbl_ltx = '\n\hline\n'.join(tbl_ltx)

os.system("echo '%s' | xclip -sel c" % tbl_ltx)
