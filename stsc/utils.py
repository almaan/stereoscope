#!/usr/bin/env python3

import torch as t
from torch.utils.data import DataLoader
import matplotlib.pyplot as plt

import numpy as np
import pandas as pd

from typing import NoReturn, List, Tuple, Union, Collection
import logging

import os.path as osp
import re
import datetime
import stsc.datasets as D
import stsc.models as M

def generate_identifier():
   return re.sub(' |:','',str(datetime.datetime.today()))


def make_joint_matrix(pths):

    index = pd.Index([])
    genes = pd.Index([])
    mlist = []
    start_pos = [0]

    for k,pth in enumerate(pths):
        cnt = read_file(pth)
        cnt.index = pd.Index([x.replace('_','&-') for x in cnt.index.values])
        mlist.append(cnt)
        index = index.append(pd.Index([str(k) + '_' + x for x in cnt.index ] ))
        genes = genes.union(cnt.columns)
        start_pos.append(cnt.shape[0])

    start_pos = np.cumsum(np.array(start_pos))
    jmat = pd.DataFrame(np.zeros((start_pos[-1],genes.shape[0])),
                       columns = genes,
                       )

    for k in range(len(start_pos) - 1):
        start = start_pos[k]
        end = start_pos[k+1] - 1
        jmat.loc[start:end,mlist[k].columns] = mlist[k].values

    jmat.columns = genes
    jmat.index = index

    return jmat

def split_joint_matrix(jmat):

    idx, name = zip(*[ idx.split('_') for idx in jmat.index ])
    name = pd.Index(name)
    idx = np.array(idx).astype(int)
    uidx = np.unique(idx)
    matlist = []

    for k in uidx:
        sel = (idx == k)
        tm = jmat.iloc[sel,:]
        tm.index = pd.Index([x.replace('&-','_') for x in name[sel].values])
        matlist.append(tm)

    return matlist

def Logger(logname,):
    """Logger for steroscope run

    Parameter
    ----------
        logname : str
            full name of file to which log should be saved

    Returns
    -------
        log : logging.Logger
            Logger object with identifier STereoSCope

    """
    log_level = logging.DEBUG

    log = logging.getLogger('stsc')
    log.setLevel(log_level)

    # filehandler to save log file
    fh = logging.FileHandler(logname)
    fh.setLevel(log_level)

    # streamhandler for stdout output
    ch = logging.StreamHandler()
    ch.setLevel(log_level)
    # format logger display message
    formatstr = '[%(asctime)s - %(name)s - %(levelname)s ] >> %(message)s'
    formatter = logging.Formatter(formatstr)
    fh.setFormatter(formatter)
    ch.setFormatter(formatter)

    log.addHandler(fh)
    log.addHandler(ch)

    return log




class SimpleProgressBar:
    """
    Progress bar to display progress during estimation

    Attributes
    ----------
        max_value : int
            total number of epochs to be used
        length: int
            number of markers to use
        symbol : str
            symbol to use as indicator

    """

    def __init__(self,
                 max_value : int,
                 length : int = 20,
                 symbol : str = "=",
                 silent_mode : bool = False,
                 ):

        self.symbol = symbol
        self.mx = max_value
        self.len = length
        self.delta = self.mx / self.len
        self.ndigits = len(str(self.mx))

        print("\r\n")

        if silent_mode:
            self.call_func = self._silent
        else:
            self.call_func = self._verbose

    def _verbose(self,
                 epoch : int,
                 value : float,
                ) -> NoReturn:

        """Updates progressbar

            Parameters
            ----------
                epoch : int
                    current epoch
                value : float
                    value to display

            """

        progress = self.symbol*int((epoch / self.delta))
        print(f"\r"
              f"Epoch : {epoch +1:<{self.ndigits}}/{self.mx:<{self.ndigits}}"
              f" | LL : {value:9E}"
              f" | \x1b[1;37m["
              f" \x1b[0;36m{progress:<{self.len}}"
              f"\x1b[1;37m]"
              f" \x1b[0m",
              end="")

    def _silent(self,
                *args,
                **kwargs,
               ) -> NoReturn:
        pass



    def __call__(self,
                 epoch : int,
                 value : float,
                 ) -> NoReturn:

        self.call_func(epoch, value)

class LossTracker:
    def __init__(self,):
        self.history = []
    def __call__(self,loss):
        self.history.append(loss)
    def __len__(self,):
        return len(self.history)
    def current(self,):
        return self.history[-1]

def read_file(file_name) :
    file = pd.read_csv(file_name,
                        header = 0,
                        index_col = 0,
                        sep = '\t')
    return file

def write_file(file,opth):
    file.to_csv(opth,
                index = True,
                header = True,
                sep = '\t')


