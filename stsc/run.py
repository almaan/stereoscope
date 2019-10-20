#!/usr/bin/env python3

import os.path as osp
from os import mkdir
import sys

import torch as t

import numpy as np
import pandas as pd

from torch.cuda import is_available
from torch.utils.data import Dataset

import stsc.fit as fit
import stsc.datasets as D
import stsc.models as M
import stsc.utils as utils
import stsc.parser as parser


def run(prs,args):

    timestamp = utils.generate_identifier()

    if len(sys.argv[1::]) < 2:
        prs.print_help()
        sys.exit(-1)

    if not osp.exists(args.out_dir):
        mkdir(args.out_dir)

    log = utils.Logger(osp.join(args.out_dir,'.'.join(['stsc',timestamp,'log'])))

    args.st_cnt = (args.st_cnt if isinstance(args.st_cnt,list) else \
                    [args.st_cnt])

    input_args = dict(sc_cnt_pth = args.sc_cnt,
                      sc_lbl_pth = args.sc_labels,
                      st_cnt_pths = args.st_cnt,
                      sc_batch_size = args.sc_batch_size,
                      topn_genes = args.topn_genes,
                      gene_list_pth = args.gene_list,
                      filter_genes = args.filter_genes,
                      min_counts = args.min_counts,
                      learning_rate = args.learning_rate,
                      st_batch_size = args.st_batch_size,
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

        R, logits,sc_model = fit.fit_sc_data(**input_args)

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

        wlist,st_model = fit.fit_st_data(R = R,
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


