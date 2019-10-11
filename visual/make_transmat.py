#!/usr/bin/env python3

import numpy as np
import pandas as pd

import os.path as osp

spotfile = "/home/alma/Documents/PhD/papers/STSC/data/hippo/st/hippo_stutility_set/spotfiles/alignment_table_Hippo3.tsv"
outdir = "/home/alma/Documents/PhD/papers/STSC/rsc/hippo_data/tmats"

def lsq_solve(arr_crd,pxl_crd):
    A = np.hstack((arr_crd.reshape(-1,1),np.ones(arr_crd.shape[0]).reshape(-1,1)))
    print(f"Shape of A is {A.shape}")
    b = pxl_crd.reshape(-1,1)

    res = np.linalg.lstsq(A,b)
    print(res)

    return res[0]


spots = pd.read_csv(spotfile,sep = '\t',header = 0, index_col = None)

print(spots.head())

solved_x = lsq_solve(spots['x'].values,spots['pixel_x'].values)
solved_y = lsq_solve(spots['y'].values,spots['pixel_y'].values)
solved_x = solved_x.flatten()
solved_y = solved_y.flatten()

print('Test functionality')
print("Test x: ")
print( spots['x'].values*solved_x[0] + solved_x[1] - spots['pixel_x'] )
print( spots['y'].values*solved_y[0] + solved_y[1] - spots['pixel_y'] )

tmat = [ float(solved_x[0]), 0, 0, 0, float(solved_y[0]), 0, float(solved_x[1]), float(solved_y[1]), 1]
tmat = [str(x) for x in tmat]
tmat = ' '.join(tmat)

opth = osp.join(outdir,''.join(['_'.join(osp.basename(spotfile).split('_')[0:-1]),'.transformation_matrix.txt']))
with open(opth, 'w+') as fopen:
    fopen.writelines(tmat)
