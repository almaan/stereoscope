#!/usr/bin/env python3

import matplotlib.pyplot as plt
import os.path as osp
import numpy as np


import sys
import time
from typing import Tuple

def rolling_average(data,
                    windowsize : int,
                   )-> np.ndarray:

    tile = np.ones(windowsize) / windowsize
    smooth = np.convolve(data,
                         tile,
                         mode = 'valid',
                        )

    return smooth

def get_loss_data(loss_file : str,
                 )-> Tuple[np.ndarray]:

    if not osp.exists(loss_file):
        print(' '.join([f"ERROR : the file {loss_file}",
                         "does not exist",
                       ]
                      ),
             )
        sys.exit(-1)

    with open(loss_file,"r+") as fopen:
        loss_history = fopen.read().lstrip(',').rstrip(',')

    loss_history = np.array([float(x) for \
                             x in loss_history.split(',')])

    epoch = np.arange(1,loss_history.shape[0]+1)

    return (epoch,
            loss_history)

def progress(loss_file : str,
             windowsize : int,
             )-> None:

    if not isinstance(windowsize,int):
        windowsize = int(windowsize)

    if windowsize % 2 == 0:
        windowsize += 1

    side = int((windowsize - 1) / 2)

    fig, ax = plt.subplots(1,
                           1,
                           figsize = (8,5),
                           num = 13,
                          )
    line1, = ax.plot([],
                    [],
                    linestyle = 'dashed',
                    color = 'black',
                  )

    line2, = ax.plot([],
                     [],
                     color = 'blue',
                     alpha = 0.2,
                     linewidth = 5,
                    )

    ax.set_ylabel('Loss',
                  fontsize = 25)
    ax.set_xlabel('Epoch',
                  fontsize = 25)

    for pos in ['top','right']:
        ax.spines[pos].set_visible(False)

    keepOn = True
    while keepOn and plt.fignum_exists(13):
        try:
            xdata,ydata = get_loss_data(loss_file)
            ydata_smooth = rolling_average(ydata,
                                           windowsize = windowsize)

            xmin,xmax = xdata.min() - 1, xdata.max() + 1
            ymin,ymax = ydata.min() - 1, ydata.max() + 1

            line1.set_xdata(xdata)
            line1.set_ydata(ydata)

            line2.set_xdata(xdata[side:-side])
            line2.set_ydata(ydata_smooth)

            ax.set_xlim([xmin,xmax])
            ax.set_ylim([ymin,ymax])

            try:
                plt.pause(10)
            except:
                keepOn = False

        except KeyboardInterrupt:
            print("Closing visualization")
            plt.close()
            keepOn = False

if __name__ == '__main__':
    progress(sys.argv[1],
             sys.argv[2])
