#!/usr/bin/env python3

import argparse as arp

def make_parser():
    prs = arp.ArgumentParser()

    parser = arp.ArgumentParser()

    parser.add_argument('-scc','--sc_cnt',
                        required = False,
                        type = str,
                        help = ''.join(["path to single cell",
                                       " count file. Should be",
                                       " on format n_cells x n_genes",
                                       " use flag sct to transpose if",
                                       " if necessary"]))

    parser.add_argument('-scl','--sc_labels',
                        required = False,
                        type = str,
                        help = ''.join(["path to single cell",
                               " labels file. Should be on",
                               ]))

    parser.add_argument('-scb','--sc_batch_size',
                        required = False,
                        default = None,
                        type = int,
                        help = ''.join(["batch size for",
                               " single cell data set",
                               ]))

    parser.add_argument('-stc','--st_cnt',
                        required = False,
                        nargs = '+',
                        help = ''.join(["path to spatial",
                               " transcriptomics count file.",
                               " Shoul be on form",
                               " n_spots x n_genes"]))

    parser.add_argument('-stm','--st_model',
                        default = None,
                        required = False,
                        help = ''.join(["path to already fitted",
                                       " st model"]))

    parser.add_argument('-scm','--sc_model',
                        required = False,
                        default = None,
                        help = ''.join(["path to already fitted",
                                       " sc model"]))



    parser.add_argument('-sce','--sc_epochs',
                        required = False,
                        default = 20000,
                        type = int,
                        help = ''.join(["number of epochs",
                                " to be used in fitting",
                                " of single cell data.",
                                " Default is set to 2e4",
                                ]))


    parser.add_argument('-stb','--st_batch_size',
                        required = False,
                        default = None,
                        type = int,
                        help = ''.join(["batch size for",
                               " st data set",
                               ]))

    parser.add_argument('-scf','--sc_fit',
                        required = False,
                        default = [None,None],
                        nargs = 2,
                        help =''.join(["parameters fitted",
                               " from single cell",
                               " data. First argument",
                               " should be path to",
                               " R-matrix and second",
                               " to logit vector"])
                               )


    parser.add_argument('-ste','--st_epochs',
                        default = 20000,
                        type = int,
                        help = ''.join(["number of epochs",
                                " to be used in fitting",
                                " of spatial transcriptomics",
                                " data.",
                                " Default is set to 2e4",
                                ]))

    parser.add_argument('-o','--out_dir',
                        required = False,
                        default = '',
                        type = str,
                        help = ''.join([" full path to output",
                                        " directory. Files will",
                                        " be saved with standard ",
                                        " name and timestamp",
                                        ]))

    parser.add_argument('-shh','--silent_mode',
                        required = False,
                        default = False,
                        action = 'store_true',
                        help = ''.join(["include to silence",
                                        "output throughout",
                                        "fitting",
                                        ]))

    parser.add_argument('-n','--topn_genes',
                        required = False,
                        default = None,
                        type = int,
                        help = ''.join(["only use top n",
                                        " mose highly expressed",
                                        " genes"
                                        ]))


    parser.add_argument('-fg','--filter_genes',
                        required = False,
                        default = False,
                        action = 'store_true',
                        help = ''.join([f"Filter Ribosomal Genes",
                                        f" and MALAT1",
                                        ]))


    parser.add_argument("-lr","--learning_rate",
                        required = False,
                        default = 0.01,
                        type = float,
                        help = ''.join([f"learning rate to be",
                                        f" used."
                                        ]))


    parser.add_argument("-mg","--min_counts",
                        required = False,
                        default = 300,
                        type = float,
                        help = ''.join([f"minimum number of ",
                                        f" counts for single cells",
                                        f" to be included in",
                                        f" the analysis",
                                        ]))


    parser.add_argument("-mc","--min_cells",
                        required = False,
                        default = 0.0,
                        type = float,
                        help = ''.join([f"minimum number of ",
                                        f" cells for genes",
                                        f" to be included in",
                                        f" the analysis",
                                        ]))


    parser.add_argument('-gp','--gpu',
                        required = False,
                        default = False,
                        action = 'store_true',
                        help = ''.join(["use gpu",
                                        ]))

    return parser
