#!/usr/bin/env python3

# Script to assemble a loompy file from the
# experiment GSE114725. Two files are required
# and provides as first and second positional arguments
# (1) file containing the raw counts as downloaded from
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE114725
# and (2) a cluster-id to cell annotation mapping with
# cluster ids in the first column and cell type annotation in
# the second.


import numpy as np
import pandas as pd
import loompy as lp

import sys
import os.path as osp

mat_pth = sys.argv[1]
annot_map_pth = sys.argv[2]

print(sys.argv)

if len(sys.argv) == 4:
    opth = osp.join(sys.argv[3],"DePeer.loom")
else:
    opth = "DePeer.loom"


mat = pd.read_csv(mat_pth, sep = ',', header = 0, index_col = None)
attrs = mat.iloc[:,0:5]
mat = mat.iloc[:,5::]


amap = pd.read_csv(annot_map_pth, sep = '\t', header = 0, index_col = None,)

mapper = { cid:bid.rstrip('_') for (cid,bid) in zip(amap.values[:,0],amap.values[:,1]) }
bio_celltype = np.array([mapper[x] for x in attrs['cluster'].values])

col_attr = {'CellID' : attrs['cellid'].values,
            'Tissue':attrs['tissue'].values,
            'Patient':attrs['patient'].values,
            'Replicate':attrs['replicate'].values,
            'ClusterID' : attrs['cluster'].values,
            'BioCelltype' : bio_celltype,
           }

row_attr = {'Gene' : mat.columns.values}

lp.create(opth, mat.values.T, row_attr, col_attr)

