#!/usr/bin/env python3

import sys
from os import mkdir
import os.path as osp
from typing import NoReturn, Union, Dict

import torch as t
from torch.utils.data import DataLoader

import numpy as np
import pandas as pd


import stsc.models as M
import stsc.datasets as D
import stsc.utils as utils

def fit(model : Union[M.ScModel,M.STModel],
        dataset : D.CountData,
        loss_tracker : utils.LossTracker,
        device : t.device,
        epochs : int,
        learning_rate : float,
        batch_size : int = None,
        silent_mode : bool = False,
        **kwargs
        ) -> None:

    """Fit Model

    Generic function to fit models

    Parameter:
    ---------

    model : Union[M.ScModel,M.STModel],
        model to be fitted
    dataset : D.CountData
        CountData dataset representing either
        single cell or ST data
    loss_tracker : utils.LossTracker
        LossTracker object
    epochs : int
        number of epochs to run
    learning_rate :
        learning rate during optimization
    batch_size : int
        batch size, if none is provided
        no batching will occur. Recommended
        to use when utilizing GPU resources
    silent_mode : bool
        whether to use silent mode or not

    Returns:
    -------
    Loss progression througout
    optimization

    """

    # move model to device
    model.to(device)
    # define optimizer
    optim = t.optim.Adam(model.parameters(),
                         lr = learning_rate)
    # instatiate progressbar
    progressBar = utils.SimpleProgressBar(epochs,
                                    silent_mode = silent_mode,
                                    length = 20)

    # use full dataset if no batch size is specified
    if batch_size is None:
        batch_size = dataset.M
    else:
        batch_size = int(np.min((batch_size,dataset.M)))

    dataloader = DataLoader(dataset,
                            batch_size = batch_size,
                            shuffle = False,
                            )

    # Use try/except to catch SIGINT
    # for early interuption
    try:
        for epoch in range(epochs):
            epoch_loss = 0.0
            for batch in dataloader:

                # move batch items to device
                for k,v in batch.items():
                    batch[k] = v.to(device)

                batch['x'].requires_grad = True
                # reset gradients
                optim.zero_grad()
                # compute loss
                loss = model.forward(**batch)
                epoch_loss += loss.item()
                # compute gradients
                loss.backward()
                # update parameters based on gradients
                optim.step()

            # update progress bar
            progressBar(epoch, epoch_loss)
            # record loss progression
            loss_tracker(epoch_loss,epoch)

        # newline after complettion
        print('\n')
        # write final loss
        loss_tracker.write_history()

    except KeyboardInterrupt:
        print(' '.join(["\n\nPress Ctrl+C again",
                        "to interrupt whole process",
                       ]
                      )
             )

def fit_st_data(st_data : D.CountData,
                R : pd.DataFrame,
                logits : pd.DataFrame,
                loss_tracker : utils.LossTracker,
                device : t.device,
                st_epochs : int,
                learning_rate : float,
                st_batch_size : int,
                silent_mode : bool = False,
                st_from_model : str = None,
                keep_noise : bool = False,
                **kwargs)->Dict[str,Union[pd.DataFrame,M.STModel]]:
    """Fit ST Data model

    Estimate proportion values for
    ST data

    Parameter:
    ---------

    st_data : D.CountData
        CountData object contining ST Data
    R :pd.DataFrame
        rates for each gene and celltype [n_genes x n_types]
    logits : pd.DataFrame
        logits for each gene [n_genes]
    loss_tracker : utils.LossTracker
    device : t.device
        device to which objects should be
        moved during optimization
    st_epochs : int
        number of epochs
    learning_rate : float
        learning rate during optimization
    st_batch_size : int
        batch_size, if non provided no
        batching will occur
    silent_mode : bool
        run in silent mode. Default False.
    st_from_model : str
        path to pre-fitted st-model state dict.
        Should be a '.pt' object
    keep_noise : bool
        keep dummy cell type in output

    Returns:
    -------
    Dictionary with estimated proportions
    and st model

    """

    # get intersection between ST data
    # and single cell parameters
    inter = st_data.intersect(R.index)

    if inter.shape[0] < 1:
        print("[ERROR] : No genes overlap in SC and"\
              " ST data. Exiting.",
              file = sys.stderr,
              )
        sys.exit(-1)

    R = R.loc[inter,:]
    logits = logits.loc[inter,:]


    t.manual_seed(1337)
    # generate ST model
    st_model = M.STModel(st_data.M,
                         R = R.values,
                         logits = logits.values,
                         device = device,
                         freeze_beta = kwargs.get("freeze_beta",False),
                         )
    # load st model from path if provided
    if st_from_model is not None:
        try:
            st_model.load_state_dict(t.load(st_from_model))
        except:
            print(' '.join(["Could not load state",
                            "dict from >> {st_from_model}"],
                          ),
                  file = sys.stderr,
                 )
    # estimate proportion values
    fit(dataset = st_data,
        model = st_model,
        loss_tracker = loss_tracker,
        device = device,
        epochs = st_epochs,
        learning_rate = learning_rate,
        batch_size = st_batch_size,
        silent_mode = silent_mode,
        )

    # get estimated unadjusted proportions
    W  = st_model.v.data.cpu().numpy().T
    # remove dummy cell type proportion values
    if not keep_noise:
        W = W[:,0:st_model.K]
        w_columns = R.columns
    else:
        w_columns = R.columns.append(pd.Index(["noise"]))
    # normalize to obtain adjusted proportions
    W = W / W.sum(axis = 1).reshape(-1,1)
    # generate pandas DataFrame from proportions
    W = pd.DataFrame(W,
                     index = st_data.index,
                     columns = w_columns)


    return {'proportions':W,
            'model':st_model,
           }

def fit_sc_data(sc_data : D.CountData,
                loss_tracker : utils.LossTracker,
                device : t.device,
                sc_epochs : int,
                sc_batch_size : int,
                learning_rate: float,
                silent_mode : bool = False,
                sc_from_model : str = None,
                **kwargs,
                )->Dict[str,Union[pd.DataFrame,
                                  M.ScModel]]:

    """Fit single cell data

    sc_data : D.CountData
        CountData Object containing
        single cell data
    loss_tracker : utils.LossTracker
    device : t.device
        device to which objects should be
        moved during optimization
    sc_epochs : int
        number of epochs
    learning_rate : float
        learning rate during optimization
    sc_batch_size : int
        batch_size, if non provided no
        batching will occur
    silent_mode : bool
        run in silent mode. Default False.
    sc_from_model : str,
        path to pre-fitted st-model state dict.
        Should be a '.pt' object

    Returns:
    -------
    Dictionary with estimated rates,
    logits and sc model

    """


    t.manual_seed(1337)
    # define single cell model
    sc_model = M.ScModel(n_genes = sc_data.G,
                         n_celltypes = sc_data.Z,
                         device = device)

    # load sc-model if provided
    if sc_from_model is not None and osp.exists(sc_from_model):
        sc_model.load_state_dict(t.load(sc_from_model))

    # fit single cell parameters
    fit(dataset = sc_data,
        model = sc_model,
        loss_tracker = loss_tracker,
        device = device,
        epochs = sc_epochs,
        learning_rate = learning_rate,
        batch_size = sc_batch_size,
        silent_mode = silent_mode
        )

    # retreive estimated parameter values
    logits = sc_model.o.data.cpu().numpy()
    R = sc_model.R.data.cpu().numpy()
    # get cell type names
    typenames = sc_data.unique_labels()

    # generate dataframes for parameters
    R = pd.DataFrame(R,
                     index = sc_data.genes,
                     columns = typenames,
                     )

    logits = pd.DataFrame(logits,
                          index = sc_data.genes,
                          columns = pd.Index(['logits']))

    return {'rates':R,
            'logits':logits,
            'model':sc_model,
           }

