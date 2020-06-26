#!/usr/bin/env python3

import re
import sys
from typing import List,Dict

import numpy as np
import pandas as pd


import torch as t
from torch.utils.data import Dataset

import stsc.utils as utils

class CountDataHelper(object):
    """
    Helper class for CountData class
    """

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
    """CountData Dataset class

    Class to hold count data from ST or
    Single Cell experiments

    Arguments:
    ---------

    cnt : pd.DataFrame
        Count data [n_observations x n_genes]
    lbl : pd.DataFrame
        Annotation/Label data [n_observations]

    """

    @CountDataHelper.update
    def __init__(self,
                 cnt : pd.DataFrame,
                 lbl : pd.DataFrame  = None,
                )-> None:

        self.cnt = cnt
        self.lbl = np.ones(self.cnt.shape[0]) * np.nan
        self.zidx = np.ones(self.cnt.shape[0]) * np.nan

        self.genes = self.cnt.columns
        self.index = self.cnt.index

        # if labels are provided
        if lbl is not None:
            self.lbl = lbl

            self.index = self.cnt.index.intersection(self.lbl.index)
            self.cnt = self.cnt.loc[self.index,:]
            self.lbl = self.lbl.loc[self.index].values.reshape(-1,)

            # convert labels to numeric indices
            tonumeric = { v:k for k,v in enumerate(np.unique(self.lbl)) }
            self.zidx = np.array([tonumeric[l] for l in self.lbl])

            # Sort data according to label enumeration
            # to speed up element acession
            srt = np.argsort(self.zidx)
            self.zidx = self.zidx[srt]
            self.lbl = self.lbl[srt]
            self.cnt = self.cnt.iloc[srt,:]

            self.zidx = t.LongTensor(self.zidx.flatten().astype(np.int32))

        # Convert to tensor
        self.cnt = t.tensor(self.cnt.values.astype(np.float32))
        self.libsize = self.cnt.sum(dim = 1)

    @CountDataHelper.update
    def filter_genes(self,
                     pattern : str = None,
                    )-> None:
        """ Filter genes based on regex-pattern

        Parameter:
        ---------
        pattern : str
            string containing regex pattern
            for which matching genes should
            be removed. Default is RP-genes
            and MALAT1

        """
        if pattern is None:
           pattern =  '^RP|MALAT1'

        keep = [ re.search(pattern,x.upper()) is \
                None for x in self.genes]

        self.cnt = self.cnt[:,keep]
        self.genes = self.genes[keep]

    @CountDataHelper.update
    def filter_bad(self,
                   min_counts : int = 0,
                   min_occurance : int = 0,
                  )-> None:

        """Filter bad data points

        Parameter:
        ----------
        min_counts : int
            minimal number of observed
            counts assigned to a specific
            spot/cell for it to be included

        min_occurance : int
            minimal number of occurances
            of a gene among all spots/cells
            for it to be included

        """

        row_thrs, col_thrs = min_counts,min_occurance
        ridx = np.where(self.cnt.sum(dim = 1) > row_thrs)[0]
        cidx = np.where((self.cnt != 0).type(t.float32).sum(dim = 0) > col_thrs)[0]

        self.cnt = self.cnt[ridx,:][:,cidx]
        self.lbl = self.lbl[ridx]
        self.zidx = self.zidx[ridx].type(t.LongTensor)

    @CountDataHelper.update
    def intersect(self,
                  exog_genes : pd.Index,
                 ) -> pd.Index:
        """Intersect genes of CountData object with external set

        Parameter:
        ----------
        exog_genes : pd.Index
            set of genes which object's gene
            set should be intersected with
        Returns:
        -------

        Pandas Index object with the intersection between
        the two sets

        """

        inter = exog_genes.intersection(self.genes)
        keep = np.array([ self.genes.get_loc(x) for x in inter])
        self.genes = inter
        self.cnt = self.cnt[:,keep]

        return self.genes

    def unique_labels(self,
                     )->np.ndarray:
        """Get unique labels

        Returns:
        -------
        Array of unique cell type labels

        """
        _,upos = np.unique(self.zidx, return_index = True)
        typenames = self.lbl[upos]

        return typenames



    def __getitem__(self,
                    idx: List[int],
                   )-> Dict:
        """Get sample with specified index

        Parameter:
        ---------
        idx : List[int]
            list of indices for samples to
            be returned

        Returns:
        -------
        Dictionary with sampe expression (x),
        label (meta), size factor (sf) and specified
        indices (gidx)

        """
        sample = {'x' : self.cnt[idx,:],
                  'meta' : self.zidx[idx],
                  'sf' : self.libsize[idx],
                  'gidx' : t.tensor(idx),
                 }

        return sample

    def __len__(self,
               )-> int:
        """Length of CountData object"""

        return self.M

def make_sc_dataset(cnt_pth : str,
                    lbl_pth : str,
                    topn_genes : int = None,
                    gene_list_pth : str = None,
                    filter_genes : bool = False,
                    lbl_colname : str = 'bio_celltype',
                    min_counts : int = 300,
                    min_cells : int = 0,
                    transpose : bool = False,
                    ):

    """
    Generate CountData object for SC-data

    Parameter:
    ---------

    cnt_pth : str
        path to SC count data

    lbl_pth : str
        path to SC label data

    topn_genes : bool
        number of top expressed genes to
        include

    gene_list_pth : str
        gene list

    lbl_colname : str
        name of column containing labels

    min_counts : int
        minimal number of observed
        counts assigned to a specific
        spot/cell for it to be included

    min_cells : int
        minimal number of occurances
        of a gene among all cells
        for it to be included

    transpose : bool
        transpose data

    Returns:
    -------

    CountData object for the SC data

    """

    sc_ext = utils.get_extenstion(cnt_pth)

    if sc_ext == 'h5ad' :
        cnt,lbl = utils.read_h5ad_sc(cnt_pth,
                                     lbl_colname,
                                     lbl_pth,
                                     )
    else:
        cnt = utils.read_file(cnt_pth,sc_ext)
        if transpose:
            cnt = cnt.T
        lbl = utils.read_file(lbl_pth)

        # get labels
        if lbl_colname is None:
            lbl = lbl.iloc[:,0]
        else:
            lbl = lbl.loc[:,lbl_colname]

    # match count and label data
    inter = cnt.index.intersection(lbl.index)
    if inter.shape[0] < 1:
        print("[ERROR] : single cell count and annotation"\
              " data did not match. Exiting.",
              file = sys.stderr,
              )
    cnt = cnt.loc[inter,:]
    lbl = lbl.loc[inter]



    # select top N expressed genes
    if topn_genes is not None:
        genesize = cnt.values.sum(axis = 0)
        topn_genes = np.min((topn_genes,genesize.shape[0]))
        sel = np.argsort(genesize)[::-1]
        sel = sel[0:topn_genes]
        cnt = cnt.iloc[:,sel]

    # only use genes in specific genes list
    # if specified
    if gene_list_pth is not None:
        with open(gene_list_pth,'r+') as fopen:
            gene_list = fopen.readlines()

        gene_list = pd.Index([ x.replace('\n','') for x in gene_list ])
        sel = cnt.columns.intersection(gene_list)
        cnt = cnt.loc[:,sel]

    # create sc data set
    dataset = CountData(cnt = cnt,
                        lbl = lbl)

    # filter genes based on names
    if filter_genes:
        dataset.filter_genes()

    # filter data based on quality
    if any([min_counts > 0,min_cells > 0]):
        dataset.filter_bad(min_counts = min_counts,
                           min_occurance = min_cells,
                          )

    return dataset


def make_st_dataset(cnt_pths : List[str],
                    topn_genes : bool = None,
                    min_counts : int = 0,
                    min_spots : int = 0,
                    filter_genes : bool = False,
                    transpose : bool = False,
                    )-> CountData :

    """
    Generate CountData object for ST-data

    Parameter:
    ---------

    cnt_pths : List[str]
        list of paths to ST-data

    topn_genes : bool
        number of top expressed genes to
        include in analysis

    min_counts : int
        minimal number of observed
        counts assigned to a specific
        spot/cell for it to be included

    min_occurance : int
        minimal number of occurances
        of a gene among all spots/cells
        for it to be included in analysis

    filter_genes : bool
        exclude MALAT1 and RP genes from
        analysis

    transpose : bool
        transpose data


    Returns:
    -------

    CountData object for the ST data

    """

    # create joint matrix for count data

    st_ext = utils.get_extenstion(cnt_pths[0])

    if st_ext == "h5ad":
        cnt = utils.read_h5ad_st(cnt_pths)
    else:
        cnt = utils.make_joint_matrix(cnt_pths,
                                      transpose)

    # select top N genes if specified
    if topn_genes is not None:
        genesize = cnt.values.sum(axis = 0)
        topn_genes = np.min((topn_genes,genesize.shape[0]))
        sel = np.argsort(genesize)[::-1]
        sel = sel[0:topn_genes]
        cnt = cnt.iloc[:,sel]

    dataset = CountData(cnt)

    # filter genes based on name
    if filter_genes:
        dataset.filter_genes()

    # filter data based on quality
    if any([min_counts > 0,min_spots > 0]):
        dataset.filter_bad(min_counts = min_counts,
                           min_occurance = min_spots,
                           )


    return dataset


