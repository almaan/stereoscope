#!/usr/bin/env python3

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import os.path as osp
from scipy.stats import pearsonr
import argparse as arp


def visualize_correlation(cmat,typenames):

    ntypes = cmat.shape[0]

    new_typenames = []
    for name in typenames:
        n = name.count('_')
        if n > 2:
            split = name.split('_')
            p1 = ' '.join(split[0:3])
            p2 = ' '.join(split[3::])
            new = p1 + '\n' + p2
            new_typenames.append(new)
        else:
            new_typenames.append(name.replace('_',' '))



    fig, ax = plt.subplots(1,1,
                           figsize = (50,40))

    cmap = plt.cm.seismic
    cmap.set_bad("gray",1)
    im = ax.imshow(cmat[:,:,0],
                   cmap = cmap,
                   aspect = 'equal')

    im.set_clim([-1,1])

    ax.set_xticks(np.arange(ntypes))
    ax.set_xticklabels(new_typenames,
                      rotation = 90,
                      fontsize = 45,
                      fontfamily = 'calibri')

    ax.set_yticks(np.arange(ntypes))
    ax.set_ylim([-0.5,ntypes-0.5])
    ax.set_yticklabels(new_typenames,
                       fontsize = 45,
                       fontfamily = 'calibri')
    cbar = fig.colorbar(im,)
    cbar.ax.tick_params(labelsize=45)

    return fig, ax



def main():
    prs = arp.ArgumentParser()

    prs.add_argument('-i','--input',
                     nargs = '+',
                     type = str,
                     required = True,
                     help = ' '.join(['Path(s) to proportion',
                                      'estimation files'])
                    )

    prs.add_argument('-o','--outdir',
                     type = str,
                     default = None,
                     required = False,
                     help = ' '.join(['Directory to save',
                                      'output'])
                    )

    prs.add_argument('-t','--tag',
                     type = str,
                     default = None,
                     required = False,
                     help = ' '.join(['tag for filename',
                                      ])
                    )

    prs.add_argument('-p','--pvalue_mask',
                     default = False,
                     action = 'store_true',
                    )

    prs.add_argument('-a','--alpha',
                     required = False,
                     default = 0.01,
                    )



    args = prs.parse_args()


    read_file = lambda file : pd.read_csv(file,
                                          sep = '\t',
                                          header = 0,
                                          index_col = 0)

    jmat = pd.DataFrame([])
    typenames = []

    if not isinstance(args.input,list):
        args.input = [args.input]

    if args.outdir is None:
        odir = osp.dirname(args.input[0])
    else:
        odir = args.outdir

    if args.tag is None:
        import datetime
        import re

        tag = re.sub(' |:','',str(datetime.datetime.today()))
    else:
        tag = args.tag


    for k,pth in enumerate(args.input):
        print(f'Reading section {pth} | {k+1} / {len(args.input)}')
        tmp = read_file(pth)
        jmat = pd.concat([jmat,tmp],
                           axis = 0)

        typenames.append(tmp.columns.values)

    is_consistent = all([all(typenames[0] == ii) for \
                         ii in typenames[1::]])

    if not is_consistent:
        print('ERROR : Proportion matrices do not have matched'
              ' cell types. Exiting')

        sys.exit(-1)

    ntypes = jmat.shape[1]
    typenames = jmat.columns.values
    cmat = np.zeros((ntypes,ntypes,2))

    for ii in range(ntypes):
        for jj in range(ii,ntypes):
            r,p = pearsonr(jmat.values[:,ii],jmat.values[:,jj])
            if args.pvalue_mask and p >= args.alpha:
                cmat[ii,jj,0] = cmat[jj,ii,0] = np.nan
            else:
                cmat[ii,jj,0] = cmat[jj,ii,0] = r

            cmat[ii,jj,1] = cmat[jj,ii,1] = p


    fig, ax = visualize_correlation(cmat,typenames)
    opth = osp.join(odir,'-'.join([tag,'correlation.png']))
    fig.tight_layout()
    fig.savefig(opth, dpi = 100)


if __name__ == '__main__':
    main()


