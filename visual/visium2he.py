#!/usr/bin/env python3


import pandas as pd
import numpy as np
import json
from PIL import Image
import matplotlib.pyplot as plt
from typing import Tuple

import argparse as arp
import os.path as osp

class ColorMap:
    def __init__(self,
                 N : int = 10,
                 cmap = plt.cm.brg,
                 push : int = 0):

        self.idxs = np.linspace(0,250,N).round().astype(int)
        self.cmap = cmap(self.idxs)
        self.N = N
        self.push = push

    def __call__(self,n):
        k = (n + self.push) % self.N
        return self.cmap[k,:]


def read_prop(pth : str) -> pd.DataFrame:
    return pd.read_csv(pth, sep = '\t',header = 0, index_col = 0)

def read_spot(pth : str ) -> pd.DataFrame:
    tmp = pd.read_csv(pth, header = None,index_col = 0)
    underTissue = tmp.values[:,0] == 1
    tmp = tmp.loc[underTissue,:]
    tmp.columns = pd.Index(['underTissue','xcoord','ycoord','xpix','ypix'])
    idx = tmp[['xcoord','ycoord']].values.astype(str)
    idx = pd.Index(['x'.join(idx[x,:]) for x in range(idx.shape[0])])
    tmp.index = idx
    return tmp

def match_data(df1 : pd.DataFrame,
               df2 : pd.DataFrame)-> Tuple[pd.DataFrame]:

    inter = df1.index.intersection(df2.index)
    df1 = df1.loc[inter,:]
    df2 = df2.loc[inter,:]

    return (df1,df2)

def get_crd(sdf : pd.DataFrame) -> np.ndarray:
    crd = sdf[['xpix','ypix']].values.astype(float)
    return crd

def json2dict(pth : str)->dict:
    with open(pth,'r+') as fopen:
        d = json.load(fopen)
    return d


def main():

    #prop_pth = "/home/alma/Documents/PhD/papers/STSC/res/hippo/hippo-visium-10x/visium_10x_1/V1_Adult_Mouse_Brain_filtered_feature_bc_matrix.tsv/W.2019-12-04071524.652591.tsv"
    #img_pth = "/home/alma/Documents/PhD/papers/STSC/rsc/hippo_visium/spatial/tissue_hires_image.L200.png"
    #spot_pth = "/home/alma/Documents/PhD/papers/STSC/rsc/hippo_visium/spatial/tissue_positions_list.csv"
    #json_pth = "/home/alma/Documents/PhD/papers/STSC/rsc/hippo_visium/spatial/scalefactors_json.json"
    #out_dir = '/tmp/hemap'

    prs = arp.ArgumentParser()

    prs.add_argument("-p","--proportions")
    prs.add_argument("-i","--image")
    prs.add_argument("-s","--spotfile")
    prs.add_argument("-j","--jsonfile")
    prs.add_argument("-o","--outdir")
    prs.add_argument("-cm","--colormap",
                     default = False,
                     action = 'store_true')

    args = prs.parse_args()

    prop = read_prop(args.proportions)
    spot = read_spot(args.spotfile)
    scale_dict = json2dict(args.jsonfile)

    prop,spot = match_data(prop,spot)
    crd = get_crd(spot)
    crd = crd * float(scale_dict['tissue_hires_scalef'])

    img = Image.open(args.image)
    img = np.asarray(img)
    he,wi,_ = img.shape

    sf = 0.4 * wi/6.4/100
    markersize = scale_dict['spot_diameter_fullres']*sf
    figsize = (wi / 100,he/100)

    cm = ColorMap(push = 2)
#    cm =  lambda x : plt.cm.Dark2( (x + 2) % 8)
    get_cluster_id = lambda x : int(x.split('_')[-1])

    for ct in prop.columns.values:
        print("[INFO] : Processing type {}".format(ct))
        fig, ax = plt.subplots(1,1,figsize = figsize)
        ax.imshow(img,aspect = 'equal')

        vals = prop[ct].values
        qv = np.quantile(vals,0.98)
        vals[qv < vals] = qv
        vals = vals / vals.max()

        #quantile scaling

        rgba = np.zeros((vals.shape[0],4))
        if args.colormap:
            rgba[:,:] = cm(get_cluster_id(ct))
            edgc = np.zeros(rgba.shape)
            edgc[:,3] = vals*0.8


        else:
            rgba[:,0] = 1
            edgc = 'none'

        rgba[:,3] = vals*0.95
        sc = ax.scatter(crd[:,1],
                        crd[:,0],
                        c = rgba,
                        s = markersize,
                        edgecolor = edgc,
        )

        tt = ax.set_title(ct)

        for pos in ax.spines.keys():
            ax.spines[pos].set_visible(False)

        ax.set_aspect("equal")
        ax.set_yticks([])
        ax.set_yticklabels([])

        ax.set_xticks([])
        ax.set_xticklabels([])

        fig.savefig(osp.join(args.outdir,ct + '.png'))
        plt.close("all")

if __name__ == '__main__':
    main()
