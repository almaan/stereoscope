#!/usr/bin/env python3

import sys
from os import mkdir, getcwd
import os.path as osp
import argparse as arp


import torch as t
from torch.cuda import is_available
from torch.utils.data import Dataset

import numpy as np
import pandas as pd


import stsc.fit as fit
import stsc.datasets as D
import stsc.models as M
import stsc.utils as utils
import stsc.parser as parser


def run(prs : arp.ArgumentParser,
        args : arp.Namespace,
       )-> None:

    """Run analysis

    Depending on specified arguments performs
    either single cell parameter estimation,
    ST-data proportion estimates or both.

    Parameter:
    ---------
    prs : argparse.ArgumentParser
    args : argparse.Namespace

    """

    # generate unique identifier for analysis
    timestamp = utils.generate_identifier()

    # ensure arguments are provided
    if len(sys.argv[1::]) < 2:
        prs.print_help()
        sys.exit(-1)

    # set output directory to cwd if none specified
    if args.out_dir is None:
        args.out_dir = getcwd()
    # create output directory if non-existant
    elif not osp.exists(args.out_dir):
        mkdir(args.out_dir)

    # instatiate logger
    log = utils.Logger(osp.join(args.out_dir,
                                '.'.join(['stsc',
                                          timestamp,
                                          'log'])
                               )
                      )

    # convert args to list if not
    args.st_cnt = (args.st_cnt if \
                   isinstance(args.st_cnt,list) else \
                   [args.st_cnt])

    # set device
    if args.gpu:
        device = t.device('cuda')
    else:
        device = t.device('cpu')

    device = (device if is_available() else t.device('cpu'))
    log.info("Using device {}".format(str(device)))

    # If parameters should be fitted from sc data
    if not all(args.sc_fit):
        log.info(' | '.join(["fitting sc data",
                             "count file : {}".format(args.sc_cnt),
                             "labels file : {}".format(args.sc_labels),
                            ])
                )

        # control that paths to sc data exists
        if not all([osp.exists(args.sc_cnt)]):

            log.error(' '.join(["One or more of the specified paths to",
                                "the sc data does not exist"]))
            sys.exit(-1)

        # load pre-fitted model if provided
        if args.sc_model is not None:
            log.info("loading state from provided sc_model")

        # Create data set for single cell data
        sc_data = D.make_sc_dataset(args.sc_cnt,
                                    args.sc_labels,
                                    topn_genes = args.topn_genes,
                                    gene_list_pth = args.gene_list,
                                    lbl_colname = args.label_colname,
                                    filter_genes = args.filter_genes,
                                    min_counts = args.min_sc_counts,
                                    min_cells = args.min_cells,
                                    transpose = args.sc_transpose,
                                    )

        log.info(' '.join(["SC data GENES : {} ".format(sc_data.G),
                           "SC data CELLS : {} ".format(sc_data.M),
                           "SC data TYPES : {} ".format(sc_data.Z),
                           ])
                )

        # generate LossTracker object
        oname_loss_track = osp.join(args.out_dir,
                                    '.'.join(["sc_loss",timestamp,"txt"])
                                   )

        sc_loss_tracker = utils.LossTracker(oname_loss_track,
                                            interval = 100,
                                            )
        # estimate parameters from single cell data
        sc_res  = fit.fit_sc_data(sc_data,
                                  loss_tracker = sc_loss_tracker,
                                  sc_epochs = args.sc_epochs,
                                  sc_batch_size = args.sc_batch_size,
                                  learning_rate = args.learning_rate,
                                  sc_from_model = args.sc_model,
                                  device = device,
                                 )

        R,logits,sc_model = sc_res['rates'],sc_res['logits'],sc_res['model']
        # save sc model
        oname_sc_model = osp.join(args.out_dir,
                                  '.'.join(['sc_model',timestamp,'pt']))

        t.save(sc_model.state_dict(),oname_sc_model)

        # save estimated parameters
        oname_R = osp.join(args.out_dir,
                           '.'.join(['R',timestamp,'tsv']))

        oname_logits = osp.join(args.out_dir,
                                '.'.join(['logits',timestamp,'tsv']))

        utils.write_file(R,oname_R)
        utils.write_file(logits,oname_logits)

    # Load already estimated single cell parameters
    elif args.st_cnt is not None:
        log.info(' | '.join(["load sc parameter",
                             "rates (R) : {}".format(args.sc_fit[0]),
                             "logodds (logits) : {}".format(args.sc_fit[1]),
                            ])
                )
        R = utils.read_file(args.sc_fit[0])
        logits = utils.read_file(args.sc_fit[1])

    # If ST data is provided estiamte proportions
    if args.st_cnt[0] is not None:
        # generate identifiying tag for each section
        sectiontag = list(map(lambda x: '.'.join(osp.basename(x).split('.')[0:-1]),args.st_cnt))
        log.info("fit st data section(s) : {}".format(args.st_cnt))

        # check that provided files exist
        if not all([osp.exists(x) for x in args.st_cnt]):
            log.error("Some of the provided ST-data paths does not exist")
            sys.exit(-1)

        if args.st_model is not None:
            log.info("loading state from provided st_model")

        # create data set for st data
        st_data =  D.make_st_dataset(args.st_cnt,
                                     topn_genes = args.topn_genes,
                                     min_counts = args.min_st_counts,
                                     min_spots = args.min_spots,
                                     filter_genes = args.filter_genes,
                                     transpose = args.st_transpose,
                                    )

        log.info(' '.join(["ST data GENES : {} ".format(st_data.G),
                           "ST data SPOTS : {} ".format(st_data.M),
                           ])
                )

        # generate LossTracker object
        oname_loss_track = osp.join(args.out_dir,
                                    '.'.join(["st_loss",timestamp,"txt"])
                                   )

        st_loss_tracker = utils.LossTracker(oname_loss_track,
                                            interval = 100,
                                           )

        # estimate proportions of cell types within st data
        st_res = fit.fit_st_data(st_data,
                                 R = R,
                                 logits = logits,
                                 loss_tracker = st_loss_tracker,
                                 st_epochs = args.st_epochs,
                                 st_batch_size = args.st_batch_size,
                                 learning_rate = args.learning_rate,
                                 silent_mode = args.silent_mode,
                                 st_from_model = args.st_model,
                                 device = device,
                                 keep_noise = args.keep_noise,
                                 freeze_beta = args.freeze_beta,
                                 )

        W,st_model = st_res['proportions'],st_res['model']

        # split joint matrix into multiple
        wlist = utils.split_joint_matrix(W)

        # save st model
        oname_st_model = osp.join(args.out_dir,
                                  '.'.join(['st_model',timestamp,'pt']))

        t.save(st_model.state_dict(),oname_st_model)

        # save st data proportion estimates results
        for s in range(len(wlist)):
            section_dir = osp.join(args.out_dir,sectiontag[s])
            if not osp.exists(section_dir):
                mkdir(section_dir)

            oname_W = osp.join(section_dir,'.'.join(['W',timestamp,'tsv']))
            log.info("saving proportions for section {} to {}".format(sectiontag[s],
                                                                      oname_W))
            utils.write_file(wlist[s],oname_W)
