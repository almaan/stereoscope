#!/usr/bin/env python3

import matplotlib.pyplot as plt
import os.path as osp
import numpy as np


import sys
import time
from typing import Tuple

def rolling_average(data : np.ndarray,
                    windowsize : int,
                   )-> np.ndarray:
    """Compute rolling average

    Parameter:
    ---------
    data : np.ndarray
        data for which rolling average should be
        computed
    windowsize : int
        size of window to use upon compuation
        of rolling average. Should be odd number.
    Retruns:
    -------
    The rolling avarege
    values for the data array.

    """

    tile = np.ones(windowsize) / windowsize
    smooth = np.convolve(data,
                         tile,
                         mode = 'valid',
                        )

    return smooth

def get_loss_data(loss_file : str,
                 )-> Tuple[np.ndarray]:
    """Read loss values from file

    Parameter:
    ---------
    loss_file : str
        path to loss file
    Returns:
    -------
    Tuple with epoch number as first
    item and loss value as second


    """
    # exit if  loss file does not exist
    if not osp.exists(loss_file):
        print(' '.join([f"ERROR : the file {loss_file}",
                         "does not exist",
                       ]
                      ),
             )
        sys.exit(-1)

    # read loss files
    with open(loss_file,"r+") as fopen:
        # remove initial and trailing commas
        loss_history = fopen.read().lstrip(',').rstrip(',')

    # convert loss history to array
    loss_history = np.array([float(x) for \
                             x in loss_history.split(',')])

    # generate epoch values
    epoch = np.arange(1,loss_history.shape[0]+1)

    return (epoch,
            loss_history)

def progress(loss_file : str,
             windowsize : int,
             )-> None:
    """Dynamic plot of loss history

    Parameter:
    ---------
    loss_file : str
        path to loss history file
    windowsize : int
        size of window to use upon compuation
        of rolling average. Should be odd number.

    """

    # make sure windowsize is int
    if not isinstance(windowsize,int):
        windowsize = int(windowsize)

    # if even windowsize value add one
    if windowsize % 2 == 0:
        windowsize += 1
    # length of array that is lost
    # in rolling average computation
    side = int((windowsize - 1) / 2)

    # create figure and axes
    fig, ax = plt.subplots(1,
                           1,
                           figsize = (8,5),
                           num = 13,
                          )
    # line to represent loss values
    line1, = ax.plot([],
                    [],
                    linestyle = 'dashed',
                    color = 'black',
                  )
    # line to represent rolling average
    # values
    line2, = ax.plot([],
                     [],
                     color = 'blue',
                     alpha = 0.2,
                     linewidth = 5,
                    )

    # customize plot
    ax.set_ylabel('Loss',
                  fontsize = 25)
    ax.set_xlabel('Epoch',
                  fontsize = 25)

    # remove spines
    for pos in ['top','right']:
        ax.spines[pos].set_visible(False)

    # update loss plot every 10th second
    keepOn = True
    while keepOn and plt.fignum_exists(13):
        try:
            # get loss data from file
            xdata,ydata = get_loss_data(loss_file)
            # compute rolling average
            ydata_smooth = rolling_average(ydata,
                                           windowsize = windowsize)

            # get limits for axes
            xmin,xmax = xdata.min() - 1, xdata.max() + 1
            ymin,ymax = ydata.min() - 1, ydata.max() + 1

            # update axes
            line1.set_xdata(xdata)
            line1.set_ydata(ydata)

            line2.set_xdata(xdata[side:-side])
            line2.set_ydata(ydata_smooth)

            ax.set_xlim([xmin,xmax])
            ax.set_ylim([ymin,ymax])

            # try except to catch interactive
            # closure (CTRL+W) of plot
            try:
                plt.pause(10)
            except:
                print("Closing visualization")
                keepOn = False

        except KeyboardInterrupt:
            print("Closing visualization")
            plt.close()
            keepOn = False

if __name__ == '__main__':
    progress(sys.argv[1],
             sys.argv[2])
