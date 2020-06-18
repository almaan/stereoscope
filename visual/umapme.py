#!/usr/bin/env python3

import sys
import os.path as osp

import numpy as np
import pandas as pd

import loompy as lp
import matplotlib.pyplot as plt


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

np.random.seed(1337)
lom_pth = "/home/alma/Documents/PhD/papers/STSC/data/hippo/sc/data/raw/l1_hippocampus.loom"
out_dir = "/tmp"
#lom_pth = sys.argv[1]
#out_dir = sys.argv[2]

ds = lp.connect(lom_pth)

cidx = np.array(['_'.join((x,str(y))) for x,y in zip(ds.ca["Class"].flatten(),ds.ca["Clusters"].flatten())])
unic,cnts = np.unique(cidx,return_counts = True)

sel = cnts > 25
keep = np.any((cidx.reshape(-1,1) == unic[sel]),axis = 1)
cidx = cidx[keep]
unic,upos,cnts = np.unique(cidx,return_index = True, return_counts = True)

cntr = np.array([[xcrd[cidx == x][0],ycrd[cidx == x][0]] for x in unic])
#cntr[:,1] +=  np.random.normal(0,0.5,size = (cntr.shape[0]))
cntr += np.array([-2.0,0.2]).reshape(1,2) * ( (-1) ** (np.arange(0,cntr.shape[0]) % 2)).reshape(-1,1)
kclusters = unic.shape[0]

cidx = ds.ca["Clusters"].flatten()[keep]
unic = cidx[upos]
classname = ds.ca["Class"].flatten()[keep][upos]



xcrd = ds.ca['_X'][keep]
ycrd = ds.ca['_Y'][keep]

#cntr = np.array([[xcrd[cidx == x].mean(),ycrd[cidx == x].mean()] for x in unic])
ncidx = unic / nclusters

#cm =  lambda x : plt.cm.Dark2( (x + 2) % 8)
cm = ColorMap(push = 2)
rgba = cm( cidx )
rgba[:,0:3] *= 1.2
rgba = np.clip(rgba,0,1)
cmap = cm(unic)

props = dict(boxstyle='square',
            # facecolor='white',
             alpha=0.8)

fig, ax = plt.subplots(1,1,figsize = (30,30))
edgecolor = np.zeros((xcrd.shape[0],4))
edgecolor[:,3] = 0.5

ax.scatter(xcrd,ycrd, c = rgba, edgecolor = edgecolor,s = 15,linewidth = 1)

for ii in range(cntr.shape[0]):
    props['facecolor'] = cmap[ii,:]
    if classname[ii][0].upper() != 'E':
        txt = classname[ii][0].upper() + unic[ii].astype(str)
    else:
        txt = classname[ii][0:2].lower().capitalize() + unic[ii].astype(str)
    print(txt)
#    ax.text(cntr[ii,0],cntr[ii,1],s = txt ,horizontalalignment = 'center',fontsize = 30,bbox = props)
#    ax.scatter([],[],c = cmap[ii,:], label = 'cluster' + str(unic[ii]),s = 80,edgecolor = 'black')

for pos in ax.spines.keys():
    ax.spines[pos].set_visible(False)


#ax.set_xlim([xcrd.min()-50,xcrd.max() + 2])
ax.set_aspect('equal')
ax.set_yticks([])
ax.set_yticklabels([])

ax.set_xticks([])
ax.set_xticklabels([])
#ax.legend(ncol = 2,loc = 'upper left',fontsize = 25)


fig.savefig(osp.join(out_dir,"mousebrain-clusters.png"))

