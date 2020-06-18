#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import umap


def umaprgb(vals,):
    rgb = umap.UMAP(n_components = 3).fit_transform(vals)

    mx = rgb.max(axis = 0 ).reshape(1,-1)
    mn = rgb.min(axis = 0 ).reshape(1,-1)

    rgb = (rgb - mn) / (mx - mn)

    return rgb


def clean_ax(axes):
    axes.set_xticks([])
    axes.set_yticks([])
    axes.set_xticklabels([])
    axes.set_yticklabels([])

    for sp in axes.spines.values():
        sp.set_visible(False)

    return axes


def val_viz(ax,
            crd,
            vals,
            markersize = 2,
            edgecolor = "none",
            cmap = None,
            alpha = None) :

    ax.scatter(crd[:,0],
               crd[:,1],
               c = vals,
               s = markersize,
               edgecolor = edgecolor,
               cmap = cmap,
               alpha = alpha,
               )

    ax = clean_ax(ax)

    ax.set_aspect("equal")

    return ax


def prop2rgb(pvals,
             color,
             alpha_max = None,
             ):

    if isinstance(color,str):
        if color[0] == '#':
            color = color.lstrip("#")
        color = list(int(color[i:i+2], 16) for i in (0, 2, 4))
        color = np.array(color) / 255.0

    rgba = np.zeros((pvals.shape[0],4))
    rgba[:,0:3] = color

    if alpha_max is not None:
        rgba[:,3] = pvals * alpha_max
    else:
        rgba[:,3] = pvals 

    return rgba

