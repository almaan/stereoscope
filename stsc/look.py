#!/usr/bin/env python3

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt

plt.rcParams.update({
    "figure.max_open_warning" : 200,
    "font.size" : 15,
    "font.family": "calibri",  # use serif/main font for text elements
})


import sys
import re
import os
import os.path as osp
import argparse as arp
import warnings
from scipy import interpolate

from sklearn.manifold import TSNE
from sklearn.decomposition import PCA

from PIL.JpegImagePlugin import JpegImageFile
import PIL.Image as Image


import umap
from numba.errors import NumbaDeprecationWarning, NumbaWarning

warnings.simplefilter('ignore', category=NumbaDeprecationWarning)
warnings.simplefilter('ignore', category=NumbaWarning)

from stsc.utils import make_joint_matrix, split_joint_matrix

#%% Funtions -------------------------------

def spltstr(string,size = 20):
    rxseps = [' ','-','\\.','_']
    if len(string) > size:
        match = re.search('|'.join(rxseps),
                          string[size::])
        if match:
            pos = size + match.start() 
            strout = spltstr(string[0:pos]) + \
                '\n' + \
                spltstr(string[pos+1::])
            return strout
        else:
            return string 
    else:
        return string


def pd2np(func):
    """Pandas to Numpy wrapper
    Wrapper to convert dataframe parameters
    to numpy array.

    """
    def wrapper(*args,**kwargs):
        cargs = list()
        ckwargs = dict()
        for arg in args:
            if isinstance(arg,pd.DataFrame):
                cargs.append(arg.values)
            else:
                cargs.append(arg)

        if len(kwargs) > 0:
            for k,v in kwargs.items():
                if isinstance(v,pd.DataFrame):
                    ckwargs.update({k:v.values})
                else:
                    ckwargs.update({k:v})

        return func(*cargs,**ckwargs)

    return wrapper

def rotation(v,theta):
    """Rotation in 3D

    Rotates a vector in 3D around axis (1,1,1)^T.
    Using Rodriguez rotatin formula :
        https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula

    Used mainly to change the hue of the compressed
    visualization.

    """
    k = np.ones((3,1)) / np.sqrt(3)
    K = np.cross(np.eye(3),k.T)
    vrot = v + (np.sin(theta)*np.dot(K,v)) + \
            (1-np.cos(theta))*np.dot(np.dot(K,K),v)

    return vrot

@pd2np
def relfreq(x, ax = 1):
    xs = x.sum(axis=ax)
    xs[xs == 0] = np.nan
    ns = ((-1,1) if ax == 1 else (1,-1))
    f =np.divide(x,xs.reshape(ns))
    f[np.isnan(f)] = 0.0
    return f


@pd2np
def ax_prop(ta1,
            x,
            y,
            pp,
            ms = 80,
            ec = "black",
            cm = plt.cm.Blues,
            mx = [35,35],
            mn = [0,0],
            alpha = 1,
            vmin = 0,
            vmax = 1):

    ta1.set_aspect('equal')
    ta1.set_xticks([])
    ta1.set_yticks([])
    ta1.set_xlim([mn[0]-1,mx[0] +1])
    ta1.set_ylim([mn[1]-1,mx[1]+1])
    ta1.autoscale(False)
    ta1.scatter(x = x,
                y = y,
                c = pp,
                cmap = cm,
                edgecolors = ec,
                s = ms,
                alpha = alpha,
                vmin = vmin,
                vmax = vmax,
                )

    for v in ta1.axes.spines.values():
                v.set_edgecolor('none')

    return ta1


def hide_spines(ax_obj):
    if not isinstance(ax_obj,np.ndarray):
        ax = np.array(ax_obj)
    else:
        ax = ax_obj.flatten()

    for p in ax:
        p.grid(False)
        p.set_xticks([])
        p.set_yticks([])
        for v in p.axes.spines.values():
                v.set_edgecolor('none')

    return ax

def compress(x,method = 'pca'):
    if method.lower() == 'tsne':
        dimred  = TSNE(n_components = 3,
                   perplexity = 20,
                   n_iter = 5000,
                   learning_rate=10,
                   n_iter_without_progress= 200,
                   )
        reduced = dimred.fit_transform(x.values)

    elif method.lower() =="umap":
        dimred = umap.UMAP(n_neighbors = 25,
                           n_components = 3)
        reduced = dimred.fit_transform(x.values)

    else:
        dimred = PCA(n_components = 3)
        reduced = dimred.fit_transform(x.values)

    reduced = rgb_transform(reduced)

    return reduced


@pd2np
def ax_compressed(ta3,x,y,v,hexagonal=False):

    xx = x.round(0).astype(int)
    yy = y.round(0).astype(int)

    minX, minY = xx.min(), yy.min()
    xx = xx - minX
    yy = yy - minY
    maxX, maxY = np.max(xx),np.max(yy)

    if not hexagonal:
        z = np.ones((maxX + 2 ,maxY + 2,3))

        for ii in range(xx.shape[0]):
            z[xx[ii]+1,yy[ii]+1] = v[ii,:]
    else :
        ncrd = np.hstack((yy.reshape(-1,1),
                          xx.reshape(-1,1)))
        z = crd2array(v,
                      ncrd,
                      maxX+2,
                      maxY+2,
                      fill = 1.0)

        mask = np.ones(z.shape)
        mask[xx,yy] = 0
        mask[xx+1,yy] = 0
        mask = mask.astype(bool)

        z[z < 0.0] = 0.0
        z[z > 1.0] = 1.0

        z[mask] = 1.0

    ta3.imshow(np.transpose(z,axes=(1,0,2)),
               interpolation = 'nearest',
               origin = 'lower',
               aspect = 'equal')

    ta3.grid(False)
    ta3.set_xticks([])
    ta3.set_yticks([])

    for v in ta3.axes.spines.values():
        v.set_edgecolor('none')

    return ta3

@pd2np
def rgb_transform(y):

    eps = 10e-12

    mn = y.min()
    mx = y.max()

    nm = (y - mn + eps ) / (mx - mn + eps)

    return nm

def read_file(pth):
    f = pd.read_csv(pth,
                    sep = '\t',
                    header = 0,
                    index_col = 0)

    return f


def get_crd(w):
    crd = [x.replace('X','').split('x') for x in w.index.values]
    crd = np.array(crd).astype(float)
    return crd

def resize_by_factor(w,h,f):
    wn = int(np.round(w/f))
    hn = int(np.round(h/f))
    return (wn,hn)

def map1d2d(s,n_cols):
    j = s % n_cols
    i = (s - j) / n_cols
    return int(i),int(j)

def crd2array(rgb,crd,w,h,fill= np.nan):
    xx = crd[:,0]
    yy = crd[:,1]
    nx = np.arange(h)
    ny = np.arange(w)
    nx, ny = np.meshgrid(nx,ny)
    arr = interpolate.griddata((xx,yy),
                               values = rgb,
                               xi = (nx,ny),
                               method = 'cubic',
                               fill_value = fill
                               )

    return arr


def look(args,):
    get_id = lambda x: '.'.join(osp.basename(x).split('.')[0:-1])
    tag = "stsc_viz"

    proppaths = args.proportions_path
    if not isinstance(proppaths,list):
        proppaths = [proppaths]

    if args.output:
        odirs = [args.output for x in range(len(proppaths))]
    else:
        odirs = [osp.join(osp.dirname(x),tag) for x in proppaths]


    basenames = [get_id(x) for x in proppaths]
    snames = [osp.basename(osp.dirname(pp)) for pp in proppaths]

    sortsynonyms = dict(section = 'section',
                        s = 'section',
                        i = "internal",
                        internal = "internal",
                        celltype = 'ct',
                        ct='ct',
                        type = 'ct')

    sort_by = sortsynonyms[args.sort_by]
    scale_by = sortsynonyms[args.scale_by]

    allwmat = make_joint_matrix(proppaths)
    allwmat[allwmat < 0] = 0.0
    allwmat.loc[:,:] = relfreq(allwmat)

    wlist = split_joint_matrix(allwmat)

    crdlist = [get_crd(w) for w in wlist]

    celltypes = allwmat.columns.tolist()
    n_sections = len(proppaths)
    n_celltypes = allwmat.shape[1]
    n_cols = args.n_cols


    # Visualize Cell Type Distribution ---------

    if sort_by == 'ct':
        n_rows = np.max((np.ceil(n_sections / n_cols).astype(int),1))
        titles = snames
        outer = n_celltypes
        inner = n_sections
        fignames = [osp.join(odirs[0],''.join([celltypes[x],'.png']))\
                    for x in range(n_celltypes)]

        suptitles = celltypes

    else:
        n_rows = np.max((np.ceil(n_celltypes / n_cols).astype(int),1))
        titles = celltypes
        outer = n_sections
        inner = n_celltypes
        fignames = [osp.join(odirs[x],''.join([snames[x],'.png'])) \
                    for x in range(n_sections)]
        suptitles = snames



    mxcrd =  [np.max(x,axis = 0) for x in crdlist]
    mncrd =  [np.min(x,axis = 0) for x in crdlist]

    figsize = ((n_cols + 1) * 3.5, (n_rows +1 ) * 3.5 + 5)

    for outside in range(outer):
        fig,ax = plt.subplots(n_rows,n_cols,figsize = figsize,squeeze = False)
        if not isinstance(ax,np.ndarray):
            ax = np.array(ax)
        else:
            ax = ax.flatten()
        for inside in range(inner):

            section_id = (inside if sort_by == 'ct' else outside )
            celltype_id = (outside if sort_by == 'ct' else inside)

            alpha = 0.01

            vmin = (np.quantile(wlist[section_id].iloc[:,celltype_id],alpha) if \
                    args.scale_by == 'i' else np.quantile(wlist[section_id].values,alpha))

            vmax = (np.quantile(wlist[section_id].iloc[:,celltype_id],1-alpha) if \
                    args.scale_by == 'i' else np.quantile(wlist[section_id].values,1-alpha))

            if args.alpha is not None:
                if args.alpha_vector:
                    alpha_vec =  wlist[section_id].iloc[:,celltype_id] * args.alpha
                else:
                    alpha_vec = args.alpha * np.ones(crdlist[section_id].shape[0])

            ax_prop(ax[inside],
                    crdlist[section_id][:,0],
                    crdlist[section_id][:,1],
                    pp = wlist[section_id].iloc[:,celltype_id],
                    mn = mncrd[section_id],
                    mx = mxcrd[section_id],
                    ms = args.marker_size,
                    ec = args.edgecolor,
                    alpha = args.alpha,
                    vmin = vmin,
                    vmax = vmax,
                    )

            ax[inside].set_title(spltstr(titles[inside]))

            if args.flip_y: ax[inside].invert_yaxis()



            hide_spines(ax[inner:ax.shape[0]])

        if not osp.exists(osp.dirname(fignames[outside])): os.mkdir(osp.dirname(fignames[outside]))
        fig.savefig(fignames[outside])

    if args.compress_method is not None:

        cmpr = compress(allwmat,
                        method = args.compress_method)

        if args.hue_rotate > 0:
            theta = np.deg2rad(args.hue_rotate)
            cmpr = rotation(cmpr.T,theta = theta).T
            cmpr = rgb_transform(cmpr)

        if args.shuffle_rgb:
            pos = np.arange(cmpr.shape[1])
            np.random.shuffle(pos)
            cmpr = cmpr[:,pos]

        cmpr = pd.DataFrame(dict(r = cmpr[:,0],
                                 g = cmpr[:,1],
                                 b = cmpr[:,2],
                                 ),
                              index = allwmat.index,
                              )

        scmpr = split_joint_matrix(cmpr)

        ct_cols = n_cols
        ct_skip = 1

        if not args.gathered_compr or len(args.proportions_path) < 2:
            for s in range(n_sections):
                figpth = osp.join(odirs[s],'.'.join([snames[s],'compressed.' + args.image_type]))
                if not osp.isdir(odirs[s]): os.mkdir(odirs[s])
                figsize = (10,10)
                fig,ax = plt.subplots(1,1,figsize = figsize)
                ax_compressed(ax,crdlist[s][:,0],crdlist[s][:,1],
                              scmpr[s],
                              hexagonal = args.hexagonal)
                if args.flip_y: ax.invert_yaxis()

                fig.tight_layout()
                fig.savefig(figpth)
                plt.close(fig)
        else:
            n_cmpr_cols = args.n_cols
            n_cmpr_rows = np.ceil(n_sections / n_cmpr_cols).astype(int)
            figpth = osp.join(odirs[0],'joint.compressed.' + args.image_type)
            if not osp.exists(odirs[0]): os.mkdir(odirs[0])

            figsize = ((n_cmpr_cols + 1) * 3.2, (n_cmpr_rows +1) * 3.2 + 1.0)

            fig, ax = plt.subplots(n_cmpr_rows,n_cmpr_cols,
                                   figsize = figsize,
                                   constrained_layout = True)

            ax = ax.reshape(n_cmpr_rows, n_cmpr_cols)

            for s in range(n_sections):
                r,c = map1d2d(s,n_cmpr_cols)
                ax_compressed(ax[r,c],
                              crdlist[s][:,0],
                              crdlist[s][:,1],
                              scmpr[s],
                              hexagonal=args.hexagonal)

                ax[r,c].set_title(spltstr(snames[s]),
                                  fontsize = 10)

                if args.flip_y: ax[r,c].invert_yaxis()

            hide_spines(ax)
            fig.savefig(figpth)

    plt.close('all')
