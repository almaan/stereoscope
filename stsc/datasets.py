#!/usr/bin/env python3

import sys
import torch as t

import numpy as np
import pandas as pd

from torch.utils.data import Dataset
from typing import List
import re

import stsc.utils as utils

class CountDataHelper(object):
    @classmethod
    def update(self,func):
        def wrapper(self,*args,**kwargs):
            tmp = func(self,*args,**kwargs)
            self.G = int(self.cnt.shape[1])
            self.M = int(self.cnt.shape[0])
            self.Z = np.unique(self.lbl).shape[0]
            self.libsize = self.cnt.sum(dim = 1)
            return tmp
        return wrapper


class CountData(Dataset):

    @CountDataHelper.update
    def __init__(self,
                 cnt : pd.DataFrame,
                 lbl : pd.DataFrame  = None,
                ):

        self.cnt = cnt
        self.lbl = np.ones(self.cnt.shape[0]) * np.nan
        self.zidx = np.ones(self.cnt.shape[0]) * np.nan

        self.genes = self.cnt.columns
        self.index = self.cnt.index

        if lbl is not None:
            self.lbl = lbl

            self.index = self.cnt.index.intersection(self.lbl.index)
            self.cnt = self.cnt.loc[self.index,:]
            self.lbl = self.lbl.loc[self.index].values.reshape(-1,)

            tonumeric = { v:k for k,v in enumerate(np.unique(self.lbl)) }
            self.zidx = np.array([tonumeric[l] for l in self.lbl])


            srt = np.argsort(self.zidx)
            self.zidx = self.zidx[srt]
            self.lbl = self.lbl[srt]
            self.cnt = self.cnt.iloc[srt,:]

            self.zidx = t.tensor(self.zidx.flatten().astype(np.int))

        self.cnt = t.tensor(self.cnt.values.astype(np.float32))
        self.libsize = self.cnt.sum(dim = 1)

    @CountDataHelper.update
    def filter_genes(self, pattern = None):
        if pattern is None:
           pattern =  '^RP|MALAT1'

        keep = [ re.search(pattern,x.upper()) is \
                None for x in self.genes]

        self.cnt = self.cnt[:,keep]
        self.genes = self.genes[keep]

    @CountDataHelper.update
    def filter_bad(self, min_counts = 300, min_cells = 0):
        row_thrs, col_thrs = min_counts,min_cells
        ridx = np.where(self.cnt.sum(dim = 1) > row_thrs)[0]
        cidx = np.where((self.cnt != 0).type(t.float32).mean(dim = 0) > col_thrs)[0]

        self.counts = self.cnt[ridx,:][:,cidx]
        self.lbl = self.lbl[ridx]
        self.zidx = self.zidx[ridx]

    @CountDataHelper.update
    def intersect(self,
                  exog_genes : pd.Index) -> pd.Index:

        inter = exog_genes.intersection(self.genes)
        keep = np.array([ self.genes.get_loc(x) for x in inter])
        self.genes = inter
        self.cnt = self.cnt[:,keep]

        return self.genes

    def unique_labels(self,):
        _,upos = np.unique(self.zidx, return_index = True)
        typenames = self.lbl[upos]
        return typenames



    def __getitem__(self, idx):
        sample = {'x' : self.cnt[idx,:],
                  'meta' : self.zidx[idx],
                  'sf' : self.libsize[idx],
                  'gidx' : t.tensor(idx),
                 }

        return sample

    def __len__(self,):
        return self.M

def make_sc_dataset(cnt_pth : str,
                    lbl_pth : str,
                    topn_genes : int = None,
                    gene_list_pth : str = None,
                    lbl_colname : str = 'bio_celltype',
                    filter_genes : bool = False,
                    min_counts : int = 0,
                    min_cells : int = 0,
                    ):


    cnt = utils.read_file(cnt_pth)
    lbl = utils.read_file(lbl_pth)

    inter = cnt.index.intersection(lbl.index)
    cnt = cnt.loc[inter,:]
    lbl = lbl.loc[inter,:]


    if lbl_colname is None:
        lbl = lbl.iloc[:,0]
    else:
        lbl = lbl.loc[:,lbl_colname]

    if topn_genes is not None:
        genesize = cnt.values.sum(axis = 0)
        topn_genes = np.min((topn_genes,genesize.shape[0]))
        sel = np.argsort(genesize)[::-1]
        sel = sel[0:topn_genes]
        cnt = cnt.iloc[:,sel]

    if gene_list is not None:
        with open(gene_list_pth,'r+') as fopen:
            gene_list = fopen.readlines()

        gene_list = pd.Index([ x.replace('\n','') for x in gene_list ])
        sel = cnt.columns.intersection(gene_list)
        cnt = cnt.iloc[:,sel]


    dataset = CountData(cnt = cnt, lbl = lbl)

    if filter_genes:
        dataset.filter_genes()

    dataset.filter_bad(min_counts, min_cells)

    return dataset


def make_st_dataset(cnt_pths : List[str],
                    topn_genes : bool = None) :

    cnt = utils.make_joint_matrix(cnt_pths)

    if topn_genes is not None:
        genesize = cnt.values.sum(axis = 0)
        topn_genes = np.min((topn_genes,genesize.shape[0]))
        sel = np.argsort(genesize)[::-1]
        sel = sel[0:topn_genes]
        cnt = cnt.iloc[:,sel]


    dataset = CountData(cnt)

    return dataset


