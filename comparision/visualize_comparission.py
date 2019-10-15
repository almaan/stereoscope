#!/usr/bin/env python3

import os
import os.path as osp
import argparse as arp

import pandas as pd
import numpy as np

import matplotlib.pyplot as plt

from typing import List,Dict

plt.rcParams.update({
        "font.family": "serif",
        "font.size": 15,
        "text.usetex": True,                    # use latex default serif font
        "font.sans-serif": ["DejaVu Sans"],  # use a specific sans-serif font
})


def make_fake_res(n_data,
                  n_types,
                  n_sets):

    def make_single(n_types,
                    n_data,
                    columns,
                    index):


        alpha = np.random.uniform(0.1,1)

        tmp = np.random.dirichlet(np.ones(n_types)*alpha, size= (n_data,))
        tmp = tmp / tmp.sum(axis = 1).reshape(-1,1)

        tmp = pd.DataFrame(tmp, index = index, columns = columns)

        return tmp

    columns = pd.Index([ "Type " + str(x) for x in range(n_types) ])
    index = pd.Index(["Spot " + str(x) for x in range(n_data) ])

    res = {'Method ' + str(x): make_single(n_types,n_data,columns,index) for \
           x in range(n_sets) }

    true = make_single(n_types,n_data,columns,index)

    return res, true

def make_values_data(n_data,
                    n_sets,
                   ):

    # to test visualization
    # and t-test

    vals = np.random.random((n_data,n_sets))
    vals = vals + np.random.uniform(0,1,size = n_sets).reshape(1,-1)
    index = pd.Index(["Spot " + str(x) for x in range(n_data) ])
    cols = pd.Index(['Method ' + str(x) for x in range(n_sets) ])
    vals = pd.DataFrame(vals, index = index, columns = cols)

    return vals

def _get_method_name(pth : str):
    return pth.split('.')[0]

def _kld_loss(estimate : np.ndarray,
              truth : np.ndarray,
              ):
    # Use estimate as Q
    # Use truth as P

    q = estimate.copy()
    p = truth.copy()

    assert np.isnan(q).sum() == 0, \
            "estimates contain zeros"

    p[p == 0] = np.nan

    terms = p*(np.log(q) - np.log(p))
    terms[np.isnan(terms) ] = 0.0
    kld = (-1.0) * terms.sum(axis = 1)

    return kld

def read_files(res_pths : List[str],
               true_pth : str,
               get_method_name = _get_method_name,
              ):
    res_list = []
    indices = []
    columns = []
    method_names = []

    for pth in res_pths:
        tmp_res = pd.read_csv(pth,
                              sep = '\t',
                              header = 0,
                              index_col = 0,
                             )

        res_list.append(tmp_res)
        indices.append(tmp_res.index.values)
        columns.append(tmp_res.columns.values)
        method_names.append(get_method_name(pth))

    indices = set.intersection(*indices)
    indices = pd.Index(list(indices))

    columns = set.intersection(*columns)
    columns = pd.Index(list(columns))


    true_prop = pd.read_csv(true_pth,
                            sep = '\t',
                            header = 0,
                            index_col = 0,
                           )

    true_prop = true_prop.loc[indices,columns]

    res_list = {method_names[x]:res_list[x].loc[indices,columns] for \
                x in range(len(res_pths)) }

    return res_list, true_prop

def generate_results_df(res_list : Dict,
                        true_prop : pd.DataFrame,
                        loss_function = _kld_loss,
                       ):

    n_spots = list(res_list.values())[0].shape[0]
    n_methods = len(res_list)

    columns = res_list.keys()
    index = list(res_list.values())[0].index

    out_res = pd.DataFrame(np.zeros((n_spots,n_methods)),
                           index = index,
                           columns = columns,
                          )

    for k,res in enumerate(res_list.values()):
        loss = loss_function(res.values,
                             true_prop.values,
                             )

        out_res.iloc[:,k] = loss

    return out_res




def visualize_results(res_data : pd.DataFrame,
                      use_raster : bool = False,
                     ):

    use_raster = True

    n_methods = res_data.shape[1]

    fig, ax = plt.subplots(1,1, figsize = (n_methods * 3,5))

    bp = ax.boxplot(res_data.T,
                    patch_artist = True,
                    zorder = 0,
                    meanline = True,
                   )


    plt.setp(bp['boxes'],
              facecolor = "#444444",
              edgecolor = 'black',
              linewidth = 2,
              )

    plt.setp(bp['medians'],
              color = 'black',
              linewidth = 2,
                )

    for pos in ['top','right']:
        ax.spines[pos].set_visible(False)


    raster_x = np.arange(1,res_data.shape[1] + 1)
    raster_x = np.random.normal(raster_x,
                                0.05,
                                size = (res_data.shape[0],raster_x.shape[0]))
    if use_raster:
        ax.scatter(raster_x,
                   res_data.values,
                   alpha = 0.1,
                   s = 10,
                   #edgecolor = 'black',
                   color = 'black',
                  )

    ax.set_xticklabels(res_data.columns,rotation = 45)
    ax.set_ylabel("Loss")
    ax.set_xlabel("Methods")

    fig.tight_layout()

    return fig, ax

if __name__ == '__main__':
    res, true = make_fake_res(n_data = 200, n_types = 5, n_sets = 8)

    sum_res = generate_results_df(res, true)

    fig, ax = visualize_results(sum_res)

    fig.savefig('/tmp/bp.png')











