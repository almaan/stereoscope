#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os.path as osp
import sys

import argparse as arp


def configure_scatter(sp,
                     size = 120.0,
                     vmin = 0,
                     vmax = 1):

    sp.set_cmap(plt.cm.Blues)
    sp.set_clim(vmin,vmax)
    sp.set_sizes(np.array([size]))
    sp.set_edgecolor((0,0,0,0.8))

def configure_ax(ax):
    ax.set_yticks([])
    ax.set_xticks([])
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.set_aspect('equal')
    for side in ax.spines.values():
        side.set_visible(False)


get_crd = lambda z: np.array([[float(x) for x in y.split('x')] for y in z])

read_file = lambda f: pd.read_csv(f,
                                  sep = '\t',
                                  index_col = 0,
                                  header = 0,
                                 )

to_list = lambda x : (x if \
                      isinstance(x,list) else \
                      [x]
                     )

prs = arp.ArgumentParser()

prs.add_argument('-w','--proportion_files',
                 type = str,
                 nargs = '+',
                 required = True,
                 help = 'Paths to proportion files',
                )

prs.add_argument('-nc','--n_columns',
                 type = int,
                 required = False,
                 default = 5,
                 help = 'number of columns',
                )

prs.add_argument('-ms','--marker_size',
                 type = int,
                 required = False,
                 default = 250,
                 help = 'size of marker in plot',
                )

prs.add_argument('-sl','--side_length',
                 type = int,
                 required = False,
                 default = 1000,
                 help = 'lenght of subplot sides',
                )
prs.add_argument('-y','--flip_y',
                 required = False,
                 default = False,
                 action = 'store_true',
                 help = 'flip y-axis',
                )

prs.add_argument('-ss','--subset_types',
                 type = str,
                 nargs = '+',
                 required = False,
                 default = None,
                 help = 'name of types to use',
                )

prs.add_argument('-o','--output_dir',
                 type = str,
                 required = False,
                 default = '.',
                 help = 'output directory',
                )

args = prs.parse_args()

print(args.output_dir)

if args.subset_types is not None:
    subset_types = to_list(args.subset_types)
else:
    subset_types = None

print(subset_types)
wpths = to_list(args.proportion_files)

if not osp.exists(args.output_dir):
    from os import mkdir
    mkdir(args.output_dir)

data = {}

for wpth in wpths:
    datum = read_file(wpth)
    if subset_types is not None:
        subset_types = datum.columns.intersection(pd.Index(subset_types))
        if subset_types.shape[0] < 1:
            print(f'ERROR: None of the specified types were present in the data',
                  f' will exit')
            sys.exit(-1)

        datum = datum.loc[:,subset_types]

    crd = get_crd(datum.index.values)
    sname = osp.basename(osp.dirname(wpth))
    print(sname)
    #vmin = np.quantile(datum.values,0.01)
    #vmax = np.quantile(datum.values,0.99)
    vmin = np.min(datum.values)
    vmax = np.max(datum.values)
    ct = datum.columns.values.tolist()
    data.update({sname:{'data':datum,
                        'coordinates':crd,
                        'vmin':vmin,
                        'vmax':vmax,
                        'celltypes':ct}
                        })

ntypes = len(ct)
wval, hval = args.side_length, args.side_length
ncols = np.min((ntypes,args.n_columns))
nrows = int(np.ceil(ntypes / ncols))

figsize = (wval/100 * ncols, hval/100 * nrows + 1.0)
marker_size = args.marker_size

for section in data.keys():
    fig, ax = plt.subplots(nrows,ncols, figsize = figsize)
    if not isinstance(ax,np.ndarray):
        ax = np.array([ax])
    ax = ax.reshape(-1,)
    for k,celltype in enumerate(data[section]['celltypes']):
        splt = ax[k].scatter(x = data[section]['coordinates'][:,0],
                             y = data[section]['coordinates'][:,1],
                             c = data[section]['data'].loc[:,celltype].values,
                             )
        ax[k].set_title(celltype)
        configure_ax(ax[k])
        configure_scatter(splt,
                          size = marker_size,
                          vmin = data[section]['vmin'],
                          vmax = data[section]['vmax'])
        if args.flip_y:
            ax[k].invert_yaxis()


    for e in range(k,ax.shape[0]):
        configure_ax(ax[e])

    oname = osp.join(args.output_dir,section + '.png')
    print(f'Writing to file : {oname}')
    fig.savefig(oname)
