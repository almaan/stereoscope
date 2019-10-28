#!/usr/bin/env python3

import os.path as osp

import pandas as pd
import numpy as np
from typing import List

import argparse as arp

def concat_files(pths : List[str],
              ):
    cnt_list = []
    indices = []
    columns = []

    for k,pth in enumerate(pths):
        tmp_res = pd.read_csv(pth,
                              sep = '\t',
                              header = 0,
                              index_col = 0,
                             )

        tmp_res.index = pd.Index(['qx' + str(k) + 'x' + x for x in tmp_res.index ])
        cnt_list.append(tmp_res)

    joint_data = pd.concat(cnt_list)

    return joint_data

def main():

    prs = arp.ArgumentParser()
    prs.add_argument('-i','--input',nargs = '+')
    prs.add_argument('-od','--outdir')
    prs.add_argument('-on','--outname')

    args = prs.parse_args()

    joint_data = concat_files(args.input)

    joint_data.to_csv(osp.join(args.outdir,args.outname + ".tsv"),
                      sep = '\t',
                      header = True,
                      index = True,
                     )

if __name__ == '__main__':
    main()


