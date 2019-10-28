#!/usr/bin/env python3

import os.path as osp
import argparse as arp

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from typing import List,Dict,NoReturn,Union


def pair_data_constant_total(estim_prop : np.ndarray,
                             truth_prop : np.ndarray,
                             truth_memb : np.ndarray,
                            ) -> Dict:

    rows,cols = np.where(truth_memb < 0)
    m_sel = truth_memb[rows,cols]


    n_memb,memb_pos = np.unique(m_sel,return_index = True)
    n_total = np.sum(np.abs(truth_memb),axis=1)[memb_pos]

    e_sel = estim_prop[rows,cols]
    t_sel = truth_prop[rows,cols]

    memb_dict = { str(np.abs(x)) + '_' + str(y):{} for \
                 x,y in zip(n_memb,n_total)}

    for key,mem in zip(memb_dict.keys(),n_memb):
        pos = (m_sel == mem)
        memb_dict[key].update({'estimate':e_sel[pos],
                               'truth':t_sel[pos]})


    return {'proportion_pairs':memb_dict,
            'n_total':n_total,
            'n_members':n_memb,
           }

def pair_data_constant_mem(estim_prop : np.ndarray,
                           truth_prop : np.ndarray,
                           truth_memb : np.ndarray,
                           ) -> Dict:

    n_total = truth_memb[0,:].sum()
    rows,cols = np.where(truth_memb < 0)
    m_sel = truth_memb[rows,cols]

    tsums = np.sum(np.abs(truth_memb),axis = 1)
    n_total,total_pos = np.unique(tsums,return_index = True)
    n_memb = np.abs(m_sel[total_pos])

    e_sel = estim_prop[rows,cols]
    t_sel = truth_prop[rows,cols]


    memb_dict = { str(np.abs(x)) + '_' + str(y):{} for x,y \
                 in zip(n_memb,n_total)}

    for key,tot in zip(memb_dict.keys(),n_total):
        pos = (tot == tsums)
        memb_dict[key].update({'estimate':e_sel[pos],
                               'truth':t_sel[pos]})


    return {'proportion_pairs':memb_dict,
            'n_total':n_total,
            'n_members':n_memb}



def hide_spines(ax,
                pos : Union[List,str] = ['top',
                                        'right'],
               ) -> NoReturn:

    if isinstance(pos,str):
        if pos == 'all':
            pos = ['top',
                   'right',
                   'left',
                   'bottom']
        else:
            pos = [pos]

    for p in pos:
        ax.spines[p].set_visible(False)

def make_hists(paired_data : Dict,
               n_cols : int = 5,
               side_length : float = 2.50,
               ):


    pairs = paired_data['proportion_pairs']
    n_total = paired_data['n_total']
    n_members = paired_data['n_members']
    n_tests = len(pairs)
    n_rows = np.ceil(n_tests / n_cols).astype(int)


    figsize = (n_cols * side_length,
               n_rows * side_length)


    fig, ax = plt.subplots(n_rows,
                           n_cols,
                           sharey = True,
                           sharex = True,
                           figsize = figsize)

    if not isinstance(ax,np.ndarray):
        ax = [ax]
    else:
        ax = ax.flatten()

    cvals = np.linspace(10,255,n_tests) / 255.0

    max_diff = {'diff':0, 'min':0,'max':0}
    for p in pairs.values():
        x = np.abs(p['estimate'] - p['truth']) / p['truth']
        x = np.sort(x)
        d = x[-1]-x[0]
        if d > max_diff['diff']:
            max_diff['diff'] = d
            max_diff['min'] = x[0]
            max_diff['max'] = x[-1]

    bins = np.linspace(max_diff['min'],
                       max_diff['max'],
                       50)

    for k,p in enumerate(pairs.values()):

        estim = p['estimate']
        truth = p['truth']
        rgb = np.zeros(3)
        rgb[0] = cvals[k]
        rgb[1] = 0.5 + ((-1)**k)*cvals[k]/2
        rgb[2] = 1 - cvals[k]

        ax[k].hist(np.abs((estim-truth)) / truth,
                   edgecolor = 'black',
                   facecolor = rgb,
                   bins = bins,
                   density = True,
                   )

        ax[k].xaxis.set_tick_params(which='both', labelbottom=True)

        ax[k].set_ylabel(r'$\frac{|w_{estim}-w_{true}|}{w_{true}}$')
       # ax[k].hist(estim,
       #            edgecolor = 'black',
       #            facecolor = rgb,
       #            bins = 50,
       #            density = True,
       #            )



        ax[k].set_title(f"{n_members[k]} / {n_total[k]} cells")
        hide_spines(ax[k])


    for ii in range(k+1,n_cols*n_rows):
        ax[ii].axis('off')


    fig.tight_layout()

    return fig,ax

def read_files(res_pths : List[str],
               true_memb_pth : str,
              ) -> Dict:
    res_list = []
    indices = []
    columns = []

    if not isinstance(res_pths,list):
        res_pths = [res_pths]

    for pth in res_pths:
        tmp_res = pd.read_csv(pth,
                              sep = '\t',
                              header = 0,
                              index_col = 0,
                             )

        res_list.append(tmp_res)
        indices.append(set(tmp_res.index.values.tolist()))
        columns.append(set(tmp_res.columns.values.tolist()))

    indices = set.intersection(*indices)
    indices = list(indices)
    indices.sort()
    indices = pd.Index(indices)

    columns = set.intersection(*columns)
    columns = list(columns)
    columns.sort()
    columns = pd.Index(columns)


    print(true_memb_pth)
    true_memb = pd.read_csv(true_memb_pth,
                            sep = '\t',
                            header = 0,
                            index_col = 0,
                           )

    true_memb = true_memb.loc[indices,columns]
    true_prop = np.abs(true_memb.values) / \
                np.abs(true_memb.values).sum(axis=1).reshape(-1,1)

    true_prop = pd.DataFrame(true_prop,
                            index = true_memb.index,
                            columns = true_memb.columns)



    res_list = [res_list[x].loc[indices,columns] for \
                x in range(len(res_pths)) ]

    return {'results':res_list,
            'truth_prop':true_prop,
            'truth_memb':true_memb,
           }


def main():
    prs = arp.ArgumentParser()

    prs.add_argument("-rf", "--result_files",
                     required = True,
                     help = "result files",
                    )

    prs.add_argument("-tf", "--truth_files",
                     required = True,
                     help = "ground truth files",
                    )

    prs.add_argument("-o","--out_dir",
                     required = True,
                     help = "",
                    )

    prs.add_argument("-c","--constant",
                     choices = ['total',
                                'member'],
                     default = 'total',
                     help = "",
                    )


    args = prs.parse_args()

    data = read_files(args.result_files,
                      args.truth_files,
                     )

    if args.constant == 'total':
        make_paired_data = pair_data_constant_total
    elif args.constant == 'member':
        make_paired_data = pair_data_constant_mem

    paired_data = make_paired_data(data['results'][0].values,
                                   data['truth_prop'].values,
                                   data['truth_memb'].values,
                                   )

    fig,ax = make_hists(paired_data)

    fig.savefig(osp.join(args.out_dir,'sensitivity.png'))


if __name__ == '__main__':
    main()


