#!/usr/bin/env python3

import os
import os.path as osp
import argparse as arp

import pandas as pd
import numpy as np
import torch as t
import torch.distributions as dists

def make_fake_dataset(n_types = 7,
                      n_cells = 40,
                      n_genes = 100):
    labels = np.repeat(np.arange(n_types).reshape(1,-1),n_cells,axis = 0).flatten()

    total_count = dists.uniform.Uniform(1,5).sample((n_types,n_genes))
    probs = dists.uniform.Uniform(0,1).sample((n_types,1))

    print(total_count)
    print(labels)

    cnt = dists.NegativeBinomial(total_count = total_count,
                                 probs = probs).sample((n_cells,))

    cnt = cnt.reshape(n_types*n_cells,n_genes)


    return cnt.numpy(), labels

def assemble_spot(cnt : np.ndarray,
                  labels : np.ndarray,
                  alpha : float = 1.0,
                  fraction : float = 0.1,
                  ):

    n_cells = dists.uniform.Uniform(low = 10,
                                    high = 30).sample().round().type(t.int)

    uni_labs, uni_counts = np.unique(labels, return_counts = True)


    assert np.all(uni_counts >=  30), \
            "Insufficient number of cells"

    n_labels = uni_labs.shape[0]

    n_types = dists.uniform.Uniform(low = 1,
                                    high =  n_labels).sample().round().type(t.int)

    pick_types = t.randperm(n_labels)[0:n_types]
    members = t.zeros(n_types)

    while members.sum() < 1:
        member_props = dists.Dirichlet(concentration = alpha * t.ones(n_types)).sample()
        members = (n_cells * member_props).round()

    props = t.zeros(n_labels)
    props[pick_types] = members / members.sum()
    members = members.type(t.int)

    spot_expr = t.zeros(cnt.shape[1]).type(t.float32)

    for z in range(n_types):
        idx = np.where(labels == uni_labs[pick_types[z]])[0]
        np.random.shuffle(idx)
        idx = idx[0:members[z]]

        # Q: should fraction be applied afterwards?
        spot_expr +=  t.tensor((cnt[idx,:]*fraction).sum(axis = 0).round().astype(np.float32))


    return spot_expr, props

def assemble_data_set(cnt : pd.DataFrame,
                      labels : pd.DataFrame,
                      n_spots : int,
                      n_genes : int,
                      **kwargs,
                     ):

    labels = labels.loc[:,'bio_celltype']

    n_genes = np.min((cnt.shape[1],n_genes))
    keep_genes = np.argsort(cnt.sum(axis=0))[::-1]
    keep_genes = keep_genes[0:n_genes]

    cnt = cnt.iloc[:,keep_genes]

    uni_labels = np.unique(labels.values)
    n_labels = uni_labels.shape[0]

    st_cnt = np.zeros((n_spots,cnt.shape[1]))
    st_prop = np.zeros((n_spots,n_labels))

    for spot in range(n_spots):
        st_cnt[spot,:], st_prop[spot,:] = assemble_spot(cnt.values,
                                                        labels.values,
                                                        **kwargs)

    index = pd.Index([ 'Spot ' + str(x + 1) for x in range(n_spots) ])

    st_cnt = pd.DataFrame(st_cnt,
                          index = index,
                          columns = cnt.columns,
                         )

    st_prop = pd.DataFrame(st_prop,
                           index = index,
                           columns = uni_labels,
                          )

    return {'counts':st_cnt,'proportions':st_prop}


def main():

    prs = arp.ArgumentParser()

    prs.add_argument('-c','--sc_counts',
                     type = str,
                    required = True)

    prs.add_argument('-l','--sc_labels',
                     type = str,
                     required = True)

    prs.add_argument('-ns','--n_spots',
                     type = int,
                     default = 1000,
                    )

    prs.add_argument('-ng','--n_genes',
                     type = int,
                     default = 500,
                    )

    prs.add_argument('-o','--out_dir',
                     default = None,
                    )

    prs.add_argument('-t','--tag',
                     default = 'st_synth',
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
