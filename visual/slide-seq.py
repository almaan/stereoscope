#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import umap

import sys

import os.path as osp
import argparse as arp

from ssplots import *


def iprint(s):
    print("[INFO] : " + s)

prs = arp.ArgumentParser()

prs.add_argument("-p","--proportions")
prs.add_argument("-o","--out_dir")
prs.add_argument("-r","--rotate",default = None, type = float)
prs.add_argument("-a","--alpha_max",default = 0.6, type = float)
prs.add_argument("-c",'--compressed',default = False, action = 'store_true')
prs.add_argument("-vp","--visualize_proportions",default = False, action = 'store_true')
prs.add_argument("-nc","--n_columns",type = int, default = 8)

args = prs.parse_args()

bname = osp.basename(args.proportions).replace('.tsv','')

iprint("Loading prop-file : {}".format(args.proportions))

props = pd.read_csv(args.proportions,
                    sep = '\t',
                    index_col = 0,
                    header = 0)

crd = np.array([x.split("x") for x in props.index]).astype(float)

if args.rotate is not None:
    iprint("Rotating {} degress".format(args.rotate))
    theta = np.deg2rad(args.rotate)
    rmat = np.array([[np.cos(theta),-np.sin(theta)],[np.sin(theta),np.cos(theta)]])
    rmat[np.abs(rmat) < 10e-10] = 0.0
    crd = np.dot(rmat,crd.T).T

types = props.columns.values.tolist()
#types = sorted(types,key = lambda x : (x.split('_')[0],
#                                       x.split('_')[1].split()[0][0],
#                                       100-len(x.split('_')[1])))
n_types = len(types)

iprint("{} types present".format(n_types))

maxpos = np.argmax(props.values,axis = 1)
maxprops = np.zeros((props.values.shape))
maxprops[np.arange(maxprops.shape[0]),maxpos] = 1
maxprops = pd.DataFrame(maxprops, columns = props.columns)


if args.visualize_proportions:
    n_cols = np.min((args.n_columns,n_types))
    n_rows =  np.ceil(n_types / n_cols).astype(int)
    figsize = (n_cols * 4.5 + 0.5, n_rows * 4.5 + 0.5)


    iprint("Assembling full proportion visualization")
    fig,ax = plt.subplots(n_rows,n_cols,
                          squeeze = False,
                          figsize = figsize,
                        )

    ax = ax.flatten()


    for k,tp in enumerate(types):

        vals = prop2rgb(maxprops[tp].values,
                        "#2E347C",
                        alpha_max = args.alpha_max,
                        )
        ax[k] = val_viz(ax[k],crd,vals,
                        markersize = 10,
                        )
        ax[k].set_title("Cluster " + tp,fontfamily = "calibri", fontsize = 35)

    for aa in ax[k+1::]:
        aa.remove()


    fig.tight_layout()
    fig.savefig(osp.join(args.out_dir,bname + "-allviz.png"))

if args.compressed:
    figsize = (20,20)
    iprint("Generating Compressed visualization")
    fig,ax = plt.subplots(1,1, figsize = figsize)
    iprint("performing umap reduction")
    np.random.seed(1337)
    vals = umaprgb(props.values)
    iprint("umap reduction complete")
    shuff = np.arange(crd.shape[0])
    np.random.shuffle(shuff)
    ax = val_viz(ax,crd[shuff,:],vals[shuff,:],markersize = 10,alpha = 0.5)
    fig.tight_layout()
    fig.savefig(osp.join(args.out_dir,bname + "-comprviz.png"))

plt.close("all")







