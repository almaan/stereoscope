#!/usr/bin/env python3

import os
import os.path as osp
import argparse as arp
from typing import List,Dict,NoReturn,Callable,Tuple
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
from rpy2.robjects import r, pandas2ri


import numpy as np
import pandas as pd
import scipy.stats as st
from scipy.stats.stats import pearsonr


import matplotlib.pyplot as plt


plt.rcParams.update({
        "font.family": "calibri",
        "font.size": 15,
})


def r_wilcox(x : np.ndarray,
             y : np.ndarray) ->  pd.DataFrame:
    """r_wilcox

    Conducts a paired one sided Wilcoxon
    signed rank test, assessing
    whether the samples in x
    in general are smaller than those
    from in y.

    Note, this function uses rbindings
    to access the slightly more
    sophisticated R implementation of 
    the test.

    Parameters
    ----------
    x : np.ndarray 
        vector of samples from population that is assumed to
        have smaller values
    y : np.ndarray 
        vector of samples from population that is assumed to
        have larger values

    Returns
    -------
    DataFrame with the results from
    the wilcoxon test.

    """

    pandas2ri.activate()

    asdf = robjects.r("as.data.frame")
    asnum = robjects.r("as.numeric")
    r_true = robjects.r("TRUE")

    a = pd.DataFrame(x)
    b = pd.DataFrame(y)
    r_a = robjects.conversion.py2rpy(a)
    r_a = asdf(r_a)
    r_b = robjects.conversion.py2rpy(b)
    r_b = asdf(r_b)
    wilcoxon = robjects.r("wilcox.test")
    args = {"y":r_b.X0,"conf.int":r_true,"alternative":"less","paired":r_true}
    pred = wilcoxon(r_a.X0,**args)
    res = dict()
    names = ['statistic','p.value','alternative','conf.int','estimate']
    for name in names:
        res.update({name:pred.rx2(name)})

    return res

def _get_method_name(pth : str,
                    ):
    return osp.basename(pth)

def _rmse(estimate : np.ndarray,
          truth : np.ndarray,
         )->np.ndarray:
    """Get RMSE

    Compute RMSE between estimate
    and true proportion values.

    RMSE is computed spotwise,
    hence one RMSE value for
    each spot its obtained

    See :
    https://en.wikipedia.org/wiki/Root-mean-square_deviation

    for information regarding RMSE

    Parameter:
    ---------
    estimate : np.ndarray
        estimated proportion values
    truth : np.ndarray

    """

    # copy values
    e = estimate.copy()
    t = truth.copy()

    # set complete zero rows to nan
    e[np.sum(e,axis=1) == 0] = np.nan
    t[np.sum(t,axis=1) == 0] = np.nan

    # normalize data
    e[np.sum(e,axis=1) == 0] = np.nan
    t[np.sum(t,axis=1) == 0] = np.nan

    # adjust nans
    e[np.isnan(e)] = 0
    t[np.isnan(t)] = 0

    # compute RMSE
    rmse = np.sqrt(((e-t)**2).mean(axis = 1))

    return rmse

def read_files(res_pths : List[str],
               true_pth : str,
               method_names :str,
              )-> Tuple[Dict[str,pd.DataFrame],pd.DataFrame]:

    """read multiple proportion files

    Parameter:
    ---------
    res_pths : List[str]
        list of paths to results files
    true_pth : str
        path to ground truth file
    method_names : str
        name of methods being compared
        should be in same order as
        res_pths

    Returns:
    -------
    Tuple where first element is
    a dictionary containing
    esimated proportion results
    and second element is
    pandas DataFrame containing
    ground truth values

    """
    res_list = []
    indices = []
    columns = []

    # read all files
    for pth in res_pths:
        tmp_res = pd.read_csv(pth,
                              sep = '\t',
                              header = 0,
                              index_col = 0,
                             )

        # keep results
        res_list.append(tmp_res)
        # build list with sets of all indices
        indices.append(set(tmp_res.index.values.tolist()))
        # build list with sets of all columns
        columns.append(set(tmp_res.columns.values.tolist()))

    # only keep spots which are found
    # in all results files
    indices = set.intersection(*indices)
    indices = list(indices)
    # sort indices
    indices.sort()
    indices = pd.Index(indices)
    # only keep types present in all
    # results files
    columns = set.intersection(*columns)
    columns = list(columns)
    columns.sort()
    columns = pd.Index(columns)

    # read ground truth files
    true_prop = pd.read_csv(true_pth,
                            sep = '\t',
                            header = 0,
                            index_col = 0,
                           )

    # only keep spots where cells are present
    is_not_zero = true_prop.index[(true_prop.values.sum(axis = 1) >  0)]
    indices = indices.intersection(is_not_zero)
    true_prop = true_prop.loc[indices,columns]
    # keep selected spots and types from results
    res_list = {method_names[x]:res_list[x].loc[indices,columns] for \
                x in range(len(res_pths)) }

    return (res_list, true_prop)

def evaluate_methods(res_list : Dict[str,pd.DataFrame],
                     true_prop : pd.DataFrame,
                     loss_function : Callable = _rmse,
                     )->pd.DataFrame:
    """Generate DataFrame of method performance

    Parameter:
    ---------
    res_list : Dict
        Dictionary of results generated by
        read_files
    true_prop : pd.DataFrame
        Data Frame with ground truth values
    loss_function :

    """

    # get number of spots present and
    # number of methods to compare
    n_spots = list(res_list.values())[0].shape[0]
    n_methods = len(res_list)

    # prepare DataFrame

    # use method names as columns
    columns = res_list.keys()
    # use spot names as indices
    index = list(res_list.values())[0].index

    out_res = pd.DataFrame(np.zeros((n_spots,
                                     n_methods)),
                           index = index,
                           columns = columns,
                          )

    for k,res in enumerate(res_list.values()):
        # compute loss for each result
        loss = loss_function(res.values,
                             true_prop.values,
                             )
        # set loss values to column representing
        # method
        out_res.loc[:,out_res.columns[k]] = loss

    return out_res


def generate_boxplot(res_data : pd.DataFrame,
                     hline : Dict[str,np.ndarray] = None,
                     use_raster : bool = False,
                     )-> Tuple[plt.Figure,
                               plt.Axes]:

    """Visualize comparison with boxplots

    Parameter:
    ---------
    res_data : pd.DataFrame
        data frame of evaluated

    Returns:
    -------

    Tuple containing Matplotlib Figure and
    Axes object

    """

    # get number of methods studied
    n_methods = res_data.shape[1]
    # generate figure and axes object
    fig, ax = plt.subplots(1,1,
                           figsize = (n_methods * 4,6))

    # if data for null distribution is
    # provided
    if hline is not None:
        ax.axhspan(hline['cil'],
                   hline['ciu'],
                   alpha = 0.2)


        ax.axhline(y = hline['mu'],
                   linestyle = 'dashed',
                   color = 'black',
                  )

    # generate boxplot
    bp = ax.boxplot(res_data.T,
                    patch_artist = True,
                    whis = 'range',
                    zorder = 0,
                   )


    plt.setp(bp['boxes'],
              facecolor = "white",
              edgecolor = 'black',
              linewidth = 1,
              )

    plt.setp(bp['medians'],
              color = 'darkred',
              linewidth = 1,
                )

    # remove top and right spines
    for pos in ['top','right']:
        ax.spines[pos].set_visible(False)


    # if raster should be plotted
    if use_raster:
        raster_x = np.arange(1,res_data.shape[1] + 1)
        raster_x = np.random.normal(raster_x,
                                0.05,
                                size = (res_data.shape[0],
                                        raster_x.shape[0]))

        ax.scatter(raster_x,
                   res_data.values,
                   alpha = 0.05,
                   s = 10,
                   color = "#444444",
                  )

    # customize graph
    ax.set_xticklabels(res_data.columns,
                       rotation = 45)

    ax.set_ylabel("RMSE")
    ax.set_xlabel("Methods")

    fig.tight_layout()

    return (fig, ax)


def generate_hist(data : pd.DataFrame,
                 )->Tuple[plt.Figure,plt.Axes]:
    """Histogram of performance metrics

    Generates histogram over the obtained
    performance metrics (eg. RMSE) for
    each method

    Parameter:
    --------
    data : pd.DataFrame
        data frame containing values
        from which histogram should
        be generated

    Returns:
    -------
    Tuple containing Matplotlib Figure and axes
    objects

    """

    fig,ax = plt.subplots(1, data.shape[1])


    if not isinstance(ax,np.ndarray):
        ax = [ax]

    for k,col in enumerate(data.columns):
        ax[k].hist(data[col].values,
                facecolor = 'gray',
                edgecolor = 'black'
               )

    return fig, ax

def test_diff(data : pd.DataFrame,
              ref_col : str,
              ) -> pd.DataFrame :
    """Conduct Wilcoxon signed-rank test

    Performs a one sided Wilcoxon signed-rank
    test; with alternative hypothesis that
    the median of the performace metric
    of competing method is less than
    that of the reference method

    Parameter:
    ---------
    data : pd.DataFrame
        pandas DataFrame generated by evaluate_methods
        function
    ref_col : str
        column name of method to be used
        as reference

    Returns:
    -------

    pandas DataFrame with average difference
    between methods and obtained p-values

    """

    # get all colnames except reference column
    non_ref = list(filter(lambda x: x != ref_col, data.columns))
    # prepare DataFrame
    columns = pd.Index(['W',
                        'estimate',
                        'p_value',
                        'CI_upper',
                        ])
    rdf = np.zeros((len(non_ref),len(columns)))
    rdf = pd.DataFrame(rdf,
                       columns = columns ,
                       index = non_ref,
                      )

    # Test reference method to others
    for k,col in enumerate(non_ref):

        # compute difference between reference and
        diff_x = data.loc[:,ref_col].values 
        diff_y = data.loc[:,col].values
        stat = r_wilcox(diff_x,diff_y)

        rdf.loc[col,'W'] = stat['statistic']
        rdf.loc[col,"CI_upper"] = stat['conf.int'][1]
        rdf.loc[col,'p_value'] = stat['p.value']
        rdf.loc[col,'estimate'] = stat['estimate']

    print(rdf)

    return rdf

def random_prop(true_prop : np.ndarray,
                niter : int = 1000,
                loss_function : Callable = _rmse,
               ) -> Dict[str,float]:
    """Random Proportion performance

    Samples proportion values from a
    Z-dimensional (Z is number of types)
    Dirichlet distribution with
    concentration parameters set to 1.0.

    Then computes the mean values
    of the performace
    metric compared to the
    ground truth niter
    times.

    Confidence interavls are taken
    as

    mean +/- 1.960*std/sqrt(niter)

    Parameter:
    ---------
    true_prop : np.ndarray
        true_proportion values
    niter : int
        number of iterantions by which
        mean performance should be computed

    loss_function : Callable
        function to compute performance metric
        from. Default is RMSE loss

    Returns:
    -------
    Dictionary with mean of performance metric,
    lower boundary and upper boundary for
    confidence interval

    """

    loss_list = []
    for iter in range(niter):
        # sample Z-dimesnional vector
        # from probability simplex
        # for every spot
        fake_prop = np.random.dirichlet(np.ones(true_prop.shape[1]),
                                        size = true_prop.shape[0])
        # compute performance metric
        loss = loss_function(fake_prop,
                             true_prop)

        # keep mean of performance metric
        loss_list.append(loss.mean())

    # compute mean and confidence
    # interval boundaries
    loss_list = np.array(loss_list)
    mu = loss_list.mean()
    std = loss_list.std()
    nsq = np.sqrt(niter)

    res = {'mu':mu,
           'cil': mu - 1.960*std/nsq,
           'ciu' : mu + 1.960*std/nsq,
          }

    return res


def main():
    """compare methods"""

    prs = arp.ArgumentParser()

    prs.add_argument("-rf", "--result_files",
                     required = True,
                     nargs = '+',
                     help = "result files",
                    )

    prs.add_argument("-tf", "--truth_file",
                     required = True,
                     help = "ground truth files",
                    )

    prs.add_argument("-o","--out_dir",
                     required = True,
                     help = "output directory",
                    )

    prs.add_argument("-mn","--method_names",
                     required = False,
                     nargs = '+',
                     default = None,
                     help = ' '.join(["name of methods",
                                      "used in comparison"]
                                    ),
                    )

    args = prs.parse_args()

    if not isinstance(args.result_files,list):
        args.result_files = [args.result_files]

    # get method names if not provided
    if args.method_names is None:
        args.method_names = [_get_method_names(x) for \
                             x in args.result_files]

    # get estimated results and true values
    res, truth = read_files(args.result_files,
                            args.truth_file,
                            args.method_names)

    # evaluate methods
    sum_res = evaluate_methods(res,
                               truth)

    # get random proportion performance
    # to use as reference
    hline = random_prop(truth.values)

    # compare methods with
    # wilcoxon test
    res_df = test_diff(sum_res, sum_res.columns[0])

    # generate boxplot of performance
    bfig, bax = generate_boxplot(sum_res,hline = hline)
    # generate histogram over loss values
    hfig,hax = generate_hist(sum_res)

    # save results
    opth_base = osp.join(args.out_dir,
                         '-'.join(args.method_names))


    bfig.savefig('-'.join([opth_base,'boxplot.png']))
    hfig.savefig('-'.join([opth_base,'histogram.png']))

    res_opth = '-'.join([opth_base,'results.txt'])

    with open(res_opth,'w+') as fopen:
        # save as latex format
        fopen.writelines(res_df.to_latex(column_format = 'c|c|c'))

if __name__ == '__main__':
    main()



