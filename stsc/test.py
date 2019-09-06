#!/usr/bin/env python3

import numpy as np

class Decorators(object):
    @classmethod
    def _update_set_info(self,func):
        def wrapper(self,*args,**kwargs):
            func(self,*args,**kwargs)
            self.G = int(self.cnt.shape[0])
            self.M = int(self.cnt.shape[1])
        return wrapper

class test:

    def __init__(self,):
        self.cnt = np.random.random((10,10))
        self.C = self.cnt[0]
        self.M = self.cnt[1]

    @Decorators._update_set_info
    def funfun(self,):
        self.cnt = self.cnt[0:4,0:6]
    def printsize(self,):
        print(self.cnt.shape)

def read_file(file_name) :
    file = pd.read_csv(counts_pth,
                        header = 0,
                        index_col = 0,
                        sep = '\t')
    return file


import pandas as pd
def make_joint_matrix(mlist):

    index = pd.Index([])
    genes = pd.Index([])
    mlist = []
    start_pos = [0]

    for k,pth in enumerate(pths):
        cnt = read_file(pth)
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
        tm.index = name[sel]
        matlist.append(tm)

    return matlist





mat1 = pd.DataFrame(np.ones((2,4)),
                    columns = ['dog','cat','peach','horse'],
                    index = ['b','c'],
                   )

mat2 = pd.DataFrame(np.ones((3,3))*4,
                    columns = ['horse','buddah','peach'],
                    index = ['a','x','c'],
                   )

Jmat = make_joint_matrix([mat1,mat2])
print('Composed')
print(Jmat)
Dmat = split_joint_matrix(Jmat)
print(Dmat[0])
print(Dmat[1])

