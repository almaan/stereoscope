#!/usr/bin/env python3

import os.path as osp
from os import mkdir
import sys

sys.path.append('/home/alma/Documents/PhD/stsc')

import torch as t

import numpy as np
import pandas as pd

from torch.cuda import is_available
from torch.utils.data import Dataset
import datasets as D
import models as M
import utils
import parser


def fit_st_data(st_cnt_pths,
                R,
                logits,
                device,
                st_epochs,
                learning_rate,
                st_batch_size,
                silent_mode,
                st_from_model,
                **kwargs):

    if not all([osp.exists(x) for x in st_cnt_pths]):
        sys.exit(-1)

    st_data = utils.make_st_dataset(st_cnt_pths)

    inter = st_data.intersect(R.index)
    R = R.loc[inter,:]
    logits = logits.loc[inter,:]

    st_model = M.STModel(st_data.M,
                         R = R.values,
                         logits = logits.values,
                         device = device)

    if st_from_model is not None and osp.exists(st_from_model):
        st_model.load_state_dict(t.load(st_from_model))

    st_loss_history = utils.fit(dataset = st_data,
                                model = st_model,
                                device = device,
                                epochs = st_epochs,
                                learning_rate = learning_rate,
                                batch_size = st_batch_size,
                                silent_mode = silent_mode,
                                )

    W  = st_model.v.data.cpu().numpy().T
    W = W[:,0:st_model.K]
    W = W / W.sum(axis = 1).reshape(-1,1)

    W = pd.DataFrame(W,
                     index = st_data.index,
                     columns = R.columns)

    wlist = utils.split_joint_matrix(W)

    return wlist,st_model

def fit_sc_data(sc_cnt_pth,
                sc_lbl_pth,
                device,
                sc_epochs,
                learning_rate,
                sc_batch_size,
                silent_mode,
                sc_from_model,
                **kwargs):

    if not osp.exists(sc_cnt_pth):
        sys.exit(-1)

    if not osp.exists(sc_lbl_pth):
        sys.exit(-1)



    sc_data = utils.make_sc_dataset(sc_cnt_pth,
                                    sc_lbl_pth,
                                    )

    sc_model = M.ScModel(n_genes = sc_data.G,
                    n_celltypes = sc_data.Z,
                    device = device)

    if sc_from_model is not None and osp.exists(sc_from_model):
        sc_model.load_state_dict(t.load(sc_from_model))



    sc_loss_history = utils.fit(dataset = sc_data,
                                model = sc_model,
                                device = device,
                                epochs = sc_epochs,
                                learning_rate = learning_rate,
                                batch_size = sc_batch_size,
                                silent_mode = silent_mode
                                )

    logits = sc_model.o.data.cpu().numpy()
    R = sc_model.R.data.cpu().numpy()

    typenames = sc_data.unique_labels()

    R = pd.DataFrame(R,
                     index = sc_data.genes,
                     columns = typenames,
                     )

    logits = pd.DataFrame(logits,
                          index = sc_data.genes,
                          columns = pd.Index(['logits']))

    return R, logits, sc_model


def main():


    timestamp = utils.generate_identifier()

    prs = parser.make_parser()
    args = prs.parse_args()

    if not osp.exists(args.out_dir):
        mkdir(args.out_dir)

    log = utils.Logger(osp.join(args.out_dir,'.'.join(['stsc',timestamp,'log'])))

    args.st_cnt = (args.st_cnt if isinstance(args.st_cnt,list) else \
                    [args.st_cnt])

    input_args = dict(sc_cnt_pth = args.sc_cnt,
                      sc_lbl_pth = args.sc_labels,
                      st_cnt_pths = args.st_cnt,
                      st_batch_size = args.sc_batch_size,
                      topn_genes = args.topn_genes,
                      filter_genes = args.filter_genes,
                      min_counts = args.min_counts,
                      learning_rate = args.learning_rate,
                      sc_batch_size = args.st_batch_size,
                      sc_epochs = args.sc_epochs,
                      st_epochs = args.st_epochs,
                      silent_mode = args.silent_mode,
                      st_from_model = args.st_model,
                      sc_from_model = args.sc_model,
                     )

    if args.gpu:
        device = t.device('cuda')
    else:
        device = t.device('cpu')

    device = (device if is_available() else t.device('cpu'))
    log.info("Using device {}".format(str(device)))
    input_args.update({'device':device})

    if not all(args.sc_fit):
        log.info("fitting sc data | count file : {} |  labels file : {}".format(args.sc_cnt,args.sc_labels))
        if args.sc_model is not None:
            log.info("loading state from provided sc_model")

        R, logits,sc_model = fit_sc_data(**input_args)

        oname_R = osp.join(args.out_dir,'.'.join(['R',timestamp,'tsv']))
        oname_logits = osp.join(args.out_dir,'.'.join(['logits',timestamp,'tsv']))
        oname_sc_model = osp.join(args.out_dir,'.'.join(['sc_model',timestamp,'pt']))

        utils.write_file(R,oname_R)
        utils.write_file(logits,oname_logits)
        t.save(sc_model.state_dict(),oname_sc_model)

    elif args.st_cnt is not None:
        log.info("load sc parameter |  rates (R) : {} | logodds (logits) : {}".format(*args.sc_fit))
        R = utils.read_file(args.sc_fit[0])
        logits = utils.read_file(args.sc_fit[1])

    if args.st_cnt is not None:

        sectiontag = list(map(lambda x: '.'.join(osp.basename(x).split('.')[0:-1]),args.st_cnt))
        log.info("fit st data section(s) : {}".format(args.st_cnt))

        if args.st_model is not None:
            log.info("loading state from provided st_model")

        wlist,st_model = fit_st_data(R = R,
                                    logits = logits,
                                     **input_args)

        oname_st_model = osp.join(args.out_dir,'.'.join(['st_model',timestamp,'pt']))
        t.save(st_model.state_dict(),oname_st_model)

        for s in range(len(wlist)):
            section_dir = osp.join(args.out_dir,sectiontag[s])
            if not osp.exists(section_dir):
                mkdir(section_dir)

            oname_W = osp.join(section_dir,'.'.join(['W',timestamp,'tsv']))
            log.info("saving proportions for section {} to {}".format(sectiontag[s],
                                                                      oname_W))
            utils.write_file(wlist[s],oname_W)


if __name__ == '__main__':
    main()



