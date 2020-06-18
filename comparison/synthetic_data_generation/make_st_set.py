#!/usr/bin/env python3

import os
import os.path as osp
import argparse as arp
from typing import Dict,Callable,List

import pandas as pd
import numpy as np
import torch as t
import torch.distributions as dists

def _assemble_spot(cnt : np.ndarray,
                  labels : np.ndarray,
                  alpha : float = 1.0,
                  fraction : float = 0.1,
                  bounds : List[int] = [10,30],
                  )->Dict[str,t.Tensor]:

    """Assemble single spot

    generates one synthetic ST-spot
    from provided single cell data

    Parameter:
    ---------
    cnt : np.ndarray
        single cell count data [n_cells x n_genes]
    labels : np.ndarray
        single cell annotations [n_cells]
    alpha : float
        dirichlet distribution
        concentration value
    fraction : float
        fraction of transcripts from each cell
        being observed in ST-spot

    Returns:
    -------
    Dictionary with expression data,
    proportion values and number of
    cells from each type at every
    spot

    """

    # sample between 10 to 30 cells to be present
    # at spot
    n_cells = dists.uniform.Uniform(low = bounds[0],
                                    high = bounds[1]).sample().round().type(t.int)

    # get unique labels found in single cell data
    uni_labs, uni_counts = np.unique(labels,
                                     return_counts = True)

    # make sure sufficient number
    # of cells are present within
    # all cell types
    assert np.all(uni_counts >=  30), \
            "Insufficient number of cells"

    # get number of different
    # cell types present
    n_labels = uni_labs.shape[0]

    # sample number of types to
    # be present at current spot
    n_types = dists.uniform.Uniform(low = 1,
                                    high =  n_labels).sample()

    n_types = n_types.round().type(t.int)

    # select which types to include
    pick_types = t.randperm(n_labels)[0:n_types]
    # pick at least one cell for spot
    members = t.zeros(n_labels).type(t.float)
    while members.sum() < 1:
        # draw proportion values from probability simplex
        member_props = dists.Dirichlet(concentration = alpha * t.ones(n_types)).sample()
        # get integer number of cells based on proportions
        members[pick_types] = (n_cells * member_props).round()

    # get proportion of each type
    props = members / members.sum()
    # convert to ints
    members = members.type(t.int)
    # get number of cells from each cell type

    # generate spot expression data
    spot_expr = t.zeros(cnt.shape[1]).type(t.float32)

    for z in range(n_types):
        # get indices of selected type
        idx = np.where(labels == uni_labs[pick_types[z]])[0]
        # pick random cells from type
        np.random.shuffle(idx)
        idx = idx[0:members[pick_types[z]]]
        # add fraction of transcripts to spot expression
        spot_expr +=  t.tensor((cnt[idx,:]*fraction).sum(axis = 0).round().astype(np.float32))


    return {'expr':spot_expr,
            'proportions':props,
            'members': members,
           }

def assemble_data_set(cnt : pd.DataFrame,
                      labels : pd.DataFrame,
                      n_spots : int,
                      n_genes : int,
                      n_cell_range : List[int],
                      assemble_fun : Callable = _assemble_spot,
                     )-> Dict[str,pd.DataFrame]:

    """Assemble Synthetic ST Data Set

    Assemble synthetic ST Data Set from
    a provided single cell data set

    Parameters:
    ----------
    cnt : pd.DataFrame
        single cell count data
    labels : pd.DataFrame
        single cell annotations
    n_spots : int
        number of spots to generate
    n_genes : int
        number of gens to include
    assemble_fun : Callable
        function to assemble single spot

    """

    # get labels
    labels = labels.loc[:,'bio_celltype']
    # make sure number of genes does not
    # exceed number of genes present
    n_genes = np.min((cnt.shape[1],n_genes))
    # select top expressed genes
    keep_genes = np.argsort(cnt.sum(axis=0))[::-1]
    keep_genes = keep_genes[0:n_genes]
    cnt = cnt.iloc[:,keep_genes]

    # get unique labels
    uni_labels = np.unique(labels.values)
    n_labels = uni_labels.shape[0]

    # prepare matrices
    st_cnt = np.zeros((n_spots,cnt.shape[1]))
    st_prop = np.zeros((n_spots,n_labels))
    st_memb = np.zeros((n_spots,n_labels))

    np.random.seed(1337)
    t.manual_seed(1337)
    # generate one spot at a time
    for spot in range(n_spots):
        spot_data = assemble_fun(cnt.values,
                                 labels.values,
                                 bounds = n_cell_range,
                                 )

        st_cnt[spot,:] = spot_data['expr']
        st_prop[spot,:] = spot_data['proportions']
        st_memb[spot,:] =  spot_data['members']

        index = pd.Index(['Spotx' + str(x + 1) for \
                          x in range(n_spots) ])

    # convert to pandas DataFrames
    st_cnt = pd.DataFrame(st_cnt,
                          index = index,
                          columns = cnt.columns,
                         )

    st_prop = pd.DataFrame(st_prop,
                           index = index,
                           columns = uni_labels,
                          )
    st_memb = pd.DataFrame(st_memb,
                           index = index,
                           columns = uni_labels,
                           )


    return {'counts':st_cnt,
            'proportions':st_prop,
            'members':st_memb}


def main():

    prs = arp.ArgumentParser()

    prs.add_argument('-c','--sc_counts',
                     type = str,
                     required = True,
                     help = ' '.join(["path to single cell",
                                      "count data",
                                     ]
                                    )
                     )

    prs.add_argument('-l','--sc_labels',
                     type = str,
                     required = True,
                     help = ' '.join(["path to single cell",
                                     "labels/annotation data",
                                     ],
                                    )

                    )

    prs.add_argument('-ns','--n_spots',
                     type = int,
                     default = 1000,
                     help = 'number of spots',
                    )

    prs.add_argument('-ng','--n_genes',
                     type = int,
                     default = 500,
                     help = 'number of genes',
                    )

    prs.add_argument('-o','--out_dir',
                     default = None,
                     help = 'output directory',
                    )

    prs.add_argument('-ncr','--n_cell_range',
                     nargs = 2,
                     default = [10,30],
                     type = int,
                     help = 'lower bound (first argument)'\
                     " and upper bound (second argument)"\
                     " for the number of cells at each spot",
                    )


    prs.add_argument('-t','--tag',
                     default = 'st_synth',
                     help = 'tag to mark data se with',
                    )

    args = prs.parse_args()

    if args.out_dir is None:
        out_dir = osp.dirname(args.sc_counts)
    else:
        out_dir = args.out_dir

    if not osp.exists(out_dir):
        os.mkdir(out_dir)

    sc_cnt_pth =  args.sc_counts
    sc_lbl_pth = args.sc_labels

    n_spots = args.n_spots
    n_genes = args.n_genes

    sc_cnt = pd.read_csv(sc_cnt_pth,
                         sep = '\t',
                         index_col = 0,
                         header = 0)

    sc_lbl = pd.read_csv(sc_lbl_pth,
                         sep = '\t',
                         index_col = 0,
                         header = 0)

    inter = sc_cnt.index.intersection(sc_lbl.index)
    sc_lbl = sc_lbl.loc[inter,:]
    sc_cnt = sc_cnt.loc[inter,:]

    assembled_set = assemble_data_set(sc_cnt,
                                      sc_lbl,
                                      n_spots = n_spots,
                                      n_genes = n_genes,
                                      n_cell_range = args.n_cell_range,
                                      assemble_fun = _assemble_spot,
                                      )

    for k,v in assembled_set.items() :
        out_pth = osp.join(out_dir, '.'.join([k,args.tag,'tsv']))
        v.to_csv(out_pth,
                    sep = '\t',
                    index = True,
                    header = True,
                   )

if __name__ == '__main__':
    main()
