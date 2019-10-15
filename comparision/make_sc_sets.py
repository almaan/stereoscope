#!/usr/bin/env python3

import os
import sys
import os.path as osp

import pandas as pd
import numpy as np
import torch as t


# Split SC data into a "validation" and "generation" set
# the generation set can be used to generate fake ST data
# whilst the validation set can be used for evaluation of the method


def main():

    sc_cnt_pth = sys.argv[1]
    sc_lbl_pth = sys.argv[2]
    out_dir = sys.argv[3]

    sc_cnt = pd.read_csv(sc_cnt_pth,
                         sep = '\t',
                         index_col = 0,
                         header = 0)

    sc_lbl = pd.read_csv(sc_lbl_pth,
                         sep = '\t',
                         index_col = 0,
                         header = 0)

    inter = sc_cnt.index.intersection(sc_lbl.index)

    sc_lbl = sc_lbl.loc[inter,:]
    sc_cnt = sc_cnt.loc[inter,:]

    labels = sc_lbl.iloc[:,0].values


    uni_labs, uni_counts = np.unique(labels,
                                     return_counts = True)

    keep_types = uni_counts > 30*2
    keep_cells = np.isin(labels, uni_labs[keep_types])

    labels = labels[keep_cells]
    sc_cnt = sc_cnt.iloc[keep_cells,:]
    sc_lbl = sc_lbl.iloc[keep_cells,:]

    uni_labs, uni_counts = np.unique(labels,
                                     return_counts = True)

    n_types = uni_labs.shape[0]

    assert np.all(uni_counts > 2), \
            "Only one cell in types"

    n_generation = (uni_counts / 2).round()
    n_validation = uni_counts - n_generation

    idx_generation = []
    idx_validation = []

    for z in range(n_types):
        tmp_idx = np.where(labels == uni_labs[z])[0]

        n_generation = int(round(tmp_idx.shape[0] / 2 ))


        idx_generation += tmp_idx[0:n_generation].tolist()
        idx_validation += tmp_idx[n_generation::].tolist()

    idx_generation.sort()
    idx_validation.sort()

    assert len(set(idx_generation).intersection(set(idx_validation))) == 0, \
            "validation and genreation set are not orthogonal"

    cnt_validation = sc_cnt.iloc[idx_validation,:]
    cnt_generation = sc_cnt.iloc[idx_generation,:]

    lbl_validation = sc_lbl.iloc[idx_validation,:]
    lbl_generation = sc_lbl.iloc[idx_generation,:]

    cnt_validation.to_csv(osp.join(out_dir,'.'.join(['validation',osp.basename(sc_cnt_pth)])),
                          sep = '\t',
                          header = True,
                          index = True,
                          index_label = 'cell',)

    cnt_generation.to_csv(osp.join(out_dir,'.'.join(['generation',osp.basename(sc_cnt_pth)])),
                          sep = '\t',
                          header = True,
                          index = True,
                          index_label = 'cell')

    lbl_validation.to_csv(osp.join(out_dir,'.'.join(['validation',osp.basename(sc_lbl_pth)])),
                          sep = '\t',
                          header = True,
                          index = True,
                          index_label = 'cell')

    lbl_generation.to_csv(osp.join(out_dir,'.'.join(['generation',osp.basename(sc_lbl_pth)])),
                          sep = '\t',
                          header = True,
                          index = True,
                          index_label = 'cell')

if __name__ == '__main__':
    main()
