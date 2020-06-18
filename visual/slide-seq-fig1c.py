#!/usr/bin/env python3

#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import umap

import sys

import os.path as osp
import argparse as arp

from ssplots import *


seldict = {
    "Neurons_59" : "#d076af",
    "Oligos_5" : "#0b73b4",
    "Neurons_22" : "#eee54c",
    "Ependymal_47" : "#59b5dc",
    "Neurons_27" : "#d8633c",
    "Vascular_68" : "#1ca174",
}

#pth = sys.argv[1]
#out_dir = sys.argv[2]

pth = "/home/alma/Documents/PhD/papers/STSC/res/hippo/slideseq_1/slide-seq-hippo.tsv/W.2019-12-10080724.760026.tsv"
out_dir = "/tmp/slide-seq-res"
use_hard = False

props = pd.read_csv(pth,
                    sep = '\t',
                    index_col = 0,
                    header = 0)

crd = np.array([x.split("x") for x in props.index]).astype(float)
rotate = 120
if use_hard:
    maxpos = np.argmax(props.values,axis = 1)
    maxprops = np.zeros((props.values.shape))
    maxprops[np.arange(maxprops.shape[0]),maxpos] = 1
    props = pd.DataFrame(maxprops, columns = props.columns, index = props.index)

if rotate is not None:
    theta = np.deg2rad(rotate)
    rmat = np.array([[np.cos(theta),-np.sin(theta)],
                     [np.sin(theta),np.cos(theta)]])
    rmat[np.abs(rmat) < 10e-10] = 0.0
    crd = np.dot(rmat,crd.T).T

figsize = (20,20)
figall, axall = plt.subplots(1,1,figsize = figsize)
for ct,col in seldict.items():
    fig, ax = plt.subplots(1,1,figsize= figsize)
    vals = prop2rgb(props[ct].values,col,alpha_max = 1)
    for xx in [ax,axall]:
        xx = val_viz(xx,crd,vals,markersize = 60,)

    ax.set_title(ct.replace('_',' '))
    fig.tight_layout()
    fig.savefig(osp.join(out_dir,ct + '.png'))

figall.tight_layout()
figall.savefig(osp.join(out_dir,"fig1c-replicate.png"))
plt.close("all")




