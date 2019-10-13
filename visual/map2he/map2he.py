#!/usr/bin/env python3

import argparse as arp
import sys
import os.path as osp
import os

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt

import PIL.Image as Image
import PIL.ImageOps

from sklearn.gaussian_process.kernels import RBF,Matern
from sklearn.gaussian_process import GaussianProcessRegressor


plt.rcParams.update({
    "font.family": "serif",  # use serif/main font for text elements
    "text.usetex": True,     # use inline math for ticks
    "pgf.rcfonts": False,    # don't setup fonts from rc parameters
    "pgf.preamble": [
         "\\usepackage{units}",          # load additional packages
         "\\usepackage{metalogo}",
         "\\usepackage{unicode-math}",   # unicode math setup
         r"\setmainfont{DejaVu Serif}",  # serif font via preamble
         ]
})

matplotlib.use('pgf')

def generate_transmat(pth : str):

    with open(pth,'r+') as fopen:
        tmat = fopen.readlines()

    tmat = [ float(x.strip('\n')) for x in tmat[0].split(' ') ]
    tmat = np.array(tmat)
    tmat = tmat.reshape(3,3)

    return tmat

def read_propmat(pth : str) :
    pmat = pd.read_csv(pth,
                       sep = '\t',
                       header = 0,
                       index_col = 0)

    return pmat

def get_crd(index : pd.Index):
    crd = [ [float(x) for x in y.replace('X','').split('x')] for y in index ]
    crd = np.array(crd)

    return crd

def get_scalemat(sf : float):
    return np.eye(2) * sf


def main():

    prs = arp.ArgumentParser()
    prs.add_argument('-i','--image',
                     required = True,
                     type = str,
                     help = 'image file',
                    )

    prs.add_argument('-t','--transformation_matrix',
                     required = True,
                     type = str,
                     help = 'transformation matrix',
                    )

    prs.add_argument('-p','--proportions_path',
                     required = True,
                     type = str,
                     help = 'output from stereoscope',
                    )

    prs.add_argument('-m','--mask_path',
                     required = False,
                     type = str,
                     default = None,
                     help = 'path to mask',
                    )


    prs.add_argument('-sf','--size_factor',
                     type = str,
                     default = '1',
                     help = ' '.join(['path to file with size factor',
                                      'or actual value'])
                    )

    prs.add_argument('-fs','--font_size',
                     type = int,
                     default = 25,
                     required = False,
                     help = ' '.join(['font size of title',])
                    )

    prs.add_argument('-ut','--use_title',
                     default = False,
                     action = 'store_true',
                     help = 'use title on images',
                    )

    prs.add_argument('-o','--odir',
                     required = False,
                     default = os.getcwd(),
                     help = 'output directory',
                    )

    prs.add_argument('-y','--flip_y',
                     required = False,
                     default = False,
                     action = 'store_true',
                     help = 'output directory',
                    )


    prs.add_argument('-g','--gray',
                     required = False,
                     default = False,
                     action = 'store_true',
                     help = 'grayscale',
                    )



    prs.add_argument('-si','--scale_internally',
                     required = False,
                     default = False,
                     action = 'store_true',
                     help = 'scale intensities internally',
                    )

    args = prs.parse_args()
    imgpth = args.image
    tmatpth = args.transformation_matrix
    wpth = args.proportions_path

    print(f"using image file {imgpth}")
    print(f"using proportions file {wpth}")

    if osp.isfile(args.size_factor):
        with open(args.size_factor,'r+') as fopen:
            resize_factor = float(fopen.readlines().strip('\n'))
    else:
        try:
            resize_factor = float(args.size_factor)
        except:
            print('[ERROR] Did not enter valid size factor path or number')
            sys.exit(-1)


    use_title = args.use_title

    tmat = generate_transmat(tmatpth)
    wmat = read_propmat(wpth)
    smat = get_scalemat(resize_factor)


    ntypes = wmat.shape[1]

    crd = np.ones((wmat.shape[0],3))
    crd[:,0:2] = get_crd(wmat.index)

    crd = np.dot(crd,tmat)[:,0:2]
    crd = np.dot(crd,smat)

    img = Image.open(imgpth)

    img_y, img_x = img.size

    if args.mask_path is not None:
        mimg = Image.open(args.mask_path).convert('RGBA')

    figsize = (img_x/ 100, img_y / 100)
    markersize = (figsize[0] * figsize[1] * 1.1 )
    adj_figsize = (figsize[0], figsize[1] + 0.1)

    id_name = osp.basename(args.image).split('_')[0]

    #color = np.array([26, 8, 82])
    #color = np.array([235, 207, 52])
    color = np.array([255,0,0])
    #color = np.array([52, 235, 79])
    color = color / color.sum()

    rgbs = np.zeros((crd.shape[0],4))
    rgbs[:,0:3] = color

    rgbe = np.zeros((crd.shape[0],4))
    rgbe[:,3] = 0.5

    for k,type in enumerate(wmat.columns.values):

        print(f"Rendering type {type}")

        dimg = np.zeros((img_y,img_x,4))


        fig,ax = plt.subplots(1,
                              1,
                              figsize = adj_figsize)
        vals = wmat[type].values

        rgbs[:,3] = vals

        if args.scale_internally:
            rgbs[:,3] /= rgbs[:,3].max()


        eimg = Image.new("RGBA",
                         (img_y,img_x),
                         color = (255,255,255,0))

        if args.mask_path is not None:
            img = Image.composite(eimg,img,mimg)

        if args.gray:
            img = img.convert('LA')

        ax.imshow(img)
        ax.scatter(crd[:,0],
                   crd[:,1],
                   marker = 'o',
                   s = markersize,
                   c = rgbs,
                   edgecolor = rgbe,
                   linewidth = 1,
                 )

        ax.set_xticks([])
        ax.set_yticks([])

        for pos in ax.spines.keys():
            ax.spines[pos].set_visible(False)

        if args.flip_y:
            ax.invert_yaxis()

        if use_title:
            y_max,y_min = crd[:,1].max(),crd[:,0].min()

            if args.flip_y:
                y_pos = y_max + 35
            else:
                y_pos = y_min - 35

            x_pos = crd[:,0].mean()

            ax.text(x_pos,
                    y_pos,
                    type,
                    verticalalignment = 'bottom',
                    horizontalalignment = 'center',
                    fontsize = args.font_size)



        opth = osp.join(args.odir,
                        ''.join([id_name,'-',type.replace(',','-'),'.png']))

        fig.tight_layout()

        fig.savefig(opth,
                    bbox_inches = 'tight',
                    pad_inches = 0,
                   transparent = True)

if __name__ == '__main__':
    main()
