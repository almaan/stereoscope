#!/usr/bin/env python3

import argparse as arp


def make_parser():
    parser = arp.ArgumentParser()

    subparsers = parser.add_subparsers(dest='command')
    run_parser = subparsers.add_parser("run",
                                       formatter_class=arp \
                                       .ArgumentDefaultsHelpFormatter)
    look_parser = subparsers.add_parser("look",
                                        formatter_class=arp \
                                        .ArgumentDefaultsHelpFormatter)
    test_parser = subparsers.add_parser("test")
    progress_parser = subparsers.add_parser('progress',
                                            formatter_class=arp \
                                            .ArgumentDefaultsHelpFormatter)

    # Run Parser Arguments ---------------------------------------------

    run_parser.add_argument('-scc', '--sc_cnt',
                            required=False,
                            type=str,
                            help='path to single cell '
                                 'count file. Should be '
                                 'in format n_cells x n_genes '
                                 'use flag sct to transpose '
                                 'if necessary')

    run_parser.add_argument('-scl', '--sc_labels',
                            required=False,
                            type=str,
                            help='path to single cell labels file. Should be on')

    run_parser.add_argument('-lcn', '--label_colname',
                            required=False,
                            default='bio_celltype',
                            type=str,
                            help='name of columns that '
                                 'cell type labels are listed')

    run_parser.add_argument('-scb', '--sc_batch_size',
                            required=False,
                            default=None,
                            type=int,
                            help='batch size for single cell data set')

    run_parser.add_argument('-stc', '--st_cnt',
                            required=False,
                            default=None,
                            nargs='+',
                            help='path to spatial '
                                 'transcriptomics count file.'
                                 'Should be in format n_spots x n_genes')

    run_parser.add_argument('-stm', '--st_model',
                            default=None,
                            required=False,
                            help='path to already fitted st model')

    run_parser.add_argument('-scm', '--sc_model',
                            required=False,
                            default=None,
                            help='path to already fitted sc model')

    run_parser.add_argument('-sce', '--sc_epochs',
                            required=False,
                            default=20000,
                            type=int,
                            help='number of epochs '
                                 'to be used in fitting '
                                 'the single cell data')

    run_parser.add_argument('-stb', '--st_batch_size',
                            required=False,
                            default=None,
                            type=int,
                            help='batch size for the spatial transcriptomics data')

    run_parser.add_argument('-scf', '--sc_fit',
                            required=False,
                            default=[None, None],
                            nargs=2,
                            help='parameters fitted '
                                 'from single cell '
                                 'data. First argument '
                                 'should be path to '
                                 'R-matrix and second '
                                 'to logit vector')

    run_parser.add_argument('-ste', '--st_epochs',
                            default=20000,
                            type=int,
                            help='number of epochs '
                                 'to be used in fitting '
                                 'the spatial transcriptomics data')

    run_parser.add_argument('-stt', '--st_transpose',
                            required=False,
                            default=False,
                            action='store_true',
                            help='transpose spatial transcriptomics data')

    run_parser.add_argument('-sct', '--sc_transpose',
                            required=False,
                            default=False,
                            action='store_true',
                            help='transpose single cell data')

    # TODO better description
    run_parser.add_argument('-kn', '--keep_noise',
                            required=False,
                            default=False,
                            action='store_true',
                            help='keep noise')

    run_parser.add_argument('-o', '--out_dir',
                            required=False,
                            default='',
                            type=str,
                            help='full path to output '
                                 'directory. Files will '
                                 'be saved with standard '
                                 'name and timestamp')

    run_parser.add_argument('-shh', '--silent_mode',
                            required=False,
                            default=False,
                            action='store_true',
                            help='include to silence '
                                 'output throughout fitting')

    run_parser.add_argument('-n', '--topn_genes',
                            required=False,
                            default=None,
                            type=int,
                            help='only use top n genes based on the '
                                 'criteria --top_criteria')

    run_parser.add_argument('-ncr', '--top_criteria',
                            required=False,
                            default='hgv',
                            choices=['hgv', 'expr'],
                            type=str,
                            help='criteria to select top genes when --topn_genes is used. '
                                 'Options are expr for expression or hgv por variance')

    run_parser.add_argument('-fg', '--filter_genes',
                            required=False,
                            default=False,
                            action='store_true',
                            help='filter ribosomal genes and Malat1')

    run_parser.add_argument('-lr', '--learning_rate',
                            required=False,
                            default=0.01,
                            type=float,
                            help='learning rate to be used in the fitting process')

    run_parser.add_argument('-mscc', '--min_sc_counts',
                            required=False,
                            default=0,
                            type=float,
                            help='minimum number of '
                                 'counts for single cells '
                                 'to be included in the analysis')

    run_parser.add_argument('-mstc', '--min_st_counts',
                            required=False,
                            default=0,
                            type=float,
                            help='minimum number of '
                                 'counts for spots '
                                 'to be included in '
                                 'the analysis')

    run_parser.add_argument('-mc', '--min_cells',
                            required=False,
                            default=0,
                            type=float,
                            help='minimum number of '
                                 'cells for genes '
                                 'to be observed in the analysis')

    run_parser.add_argument('-ms', '--min_spots',
                            required=False,
                            default=0,
                            type=float,
                            help='minimum number of '
                                 'spots for genes '
                                 'to be observed in '
                                 'the analysis')

    run_parser.add_argument('-gp', '--gpu',
                            required=False,
                            default=False,
                            action='store_true',
                            help='use gpu accelerated computation')

    run_parser.add_argument('-gl', '--gene_list',
                            required=False,
                            default=None,
                            type=str,
                            help='path to list of genes to use in the analysis')

    run_parser.add_argument('-sub', '--sc_upper_bound',
                            required=False,
                            default=None,
                            type=int,
                            help='upper bound limit for single cell subsampling')

    run_parser.add_argument('-slb', '--sc_lower_bound',
                            required=False,
                            default=None,
                            type=int,
                            help='lower bound limit for single cell subsampling')

    # TODO better description
    run_parser.add_argument('-fb', '--freeze_beta',
                            default=False,
                            action='store_true',
                            help='freeze beta parameter')

    # Look Parser Arguments -----------------------------------------------

    look_parser.add_argument('-pp', '--proportions_path',
                             type=str,
                             nargs='+',
                             required=True,
                             help='path to proportions '
                                  'file generated by '
                                  'stereoscope run. Named W*.tsv'
                             )
    # TODO add choices here
    look_parser.add_argument('-c', '--compress_method',
                             type=str,
                             required=False,
                             default=None,
                             help='method to be used '
                                  'for compression of '
                                  'information'
                             )

    look_parser.add_argument('-ms', '--marker_size',
                             type=int,
                             required=False,
                             default=100,
                             help='size of scatter plot markers'
                             )

    look_parser.add_argument('-o', '--output',
                             type=str,
                             required=False,
                             default='',
                             help='path to output '
                                  'can either be '
                                  'a directory or a '
                                  'filename. If only '
                                  'a dir is given, same '
                                  'basename the the plots is used'
                             )

    look_parser.add_argument('-nc', '--n_cols',
                             default=2,
                             type=int,
                             required=False
                             )

    look_parser.add_argument('-y', '--flip_y',
                             required=False,
                             default=False,
                             action='store_true'
                             )

    look_parser.add_argument('-sb', '--sort_by',
                             required=False,
                             default='ct',
                             type=str
                             )

    look_parser.add_argument('-sc', '--scale_by',
                             required=False,
                             default='ct',
                             type=str
                             )

    look_parser.add_argument('-gb', '--gathered_compr',
                             required=False,
                             default=False,
                             action='store_true'
                             )

    look_parser.add_argument('-sf', '--scaling_factor',
                             required=False,
                             type=float,
                             default=4.0
                             )

    look_parser.add_argument('-hu', '--hue_rotate',
                             required=False,
                             type=float,
                             default=-1.0
                             )

    look_parser.add_argument('-hex', '--hexagonal',
                             required=False,
                             default=False,
                             action='store_true'
                             )

    look_parser.add_argument('-ec', '--edgecolor',
                             required=False,
                             type=str,
                             default='black',
                             help="spot edgecolor"
                             )

    look_parser.add_argument('-ext', '--image_type',
                             required=False,
                             type=str,
                             default='png',
                             help='image file type'
                             )

    look_parser.add_argument('-al', '--alpha',
                             required=False,
                             type=float,
                             default=1,
                             help='facecolor alpha',
                             )

    look_parser.add_argument('-av', '--alpha_vector',
                             required=False,
                             default=False,
                             action='store_true',
                             help='use value based alpha',
                             )

    look_parser.add_argument('-cm', '--colormap',
                             required=False,
                             default="Blues",
                             type=str,
                             help='name of matplotlib colormap to use'
                             )

    look_parser.add_argument('-thr', '--threshold',
                             required=False,
                             default=None,
                             type=float,
                             help='threshold value for proportion visualization'
                             )

    look_parser.add_argument('-ht', '--hard_type',
                             required=False,
                             default=False,
                             action="store_true",
                             help='make hard type plot'
                             )

    look_parser.add_argument('-io', '--image_orientation',
                             required=False,
                             default=False,
                             action="store_true",
                             help='arrange capture locations '
                                  'in same orientation as HE image'
                             )

    look_parser.add_argument('-ss', '--side_size',
                             required=False,
                             default=350,
                             type=float,
                             help='subplot side size',
                             )

    look_parser.add_argument('-shu', '--shuffle_rgb',
                             required=False,
                             default=False,
                             action='store_true',
                             help='shuffle RGB colors '
                                  'in the compressed '
                                  'visualization'
                             )

    # Progress Parser Arguments -----------------------------------------------

    progress_parser.add_argument('-lf', '--loss_file',
                                 required=True,
                                 help='path to loss data'
                                 )

    progress_parser.add_argument('-ws', '--window_size',
                                 required=False,
                                 default=11,
                                 help='window size for rolling average'
                                 )

    return parser
