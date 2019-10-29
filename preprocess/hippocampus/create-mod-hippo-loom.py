#!/usr/bin/env python3

import os
import sys
import os.path as osp

import loompy as lp
import numpy as np

pth = sys.argv[1]
odir = sys.argv[2]

ds = lp.connect(pth)

celltype_1 = np.array([ '_'.join([str(x).replace(',','_'),str(y)]) for \
                       x,y in zip(ds.ca['Class'],ds.ca['Clusters']) ])

celltype_2 = np.array([ '_'.join([str(x).replace(',','_'),str(y)]) for \
                       x,y in zip(ds.ca['Subclass'],ds.ca['Clusters']) ])


new_ds_ca = dict(CellID = ds.ca['CellID'],
                 Clusters = ds.ca['Clusters'],
                 Subclass = ds.ca['Subclass'],
                 Celltype_1 = celltype_1,
                 Celltype_2 = celltype_2)

new_ds_ra = dict(Gene = ds.ra['Gene'])
new_basename = '_'.join(['mod',osp.basename(pth)])
filename = os.path.join(odir,new_basename)
lp.create(filename, ds[:,:], new_ds_ra, new_ds_ca)
print(f"successfully created modifed loom-file >> {new_basename}")
