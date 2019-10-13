#!/usr/bin/env python3

import os.path as osp
from os import mkdir
import sys

from typing import NoReturn

import torch as t
import numpy as np

from torch.utils.data import DataLoader
import stsc.utils

def fit(model,
        dataset,
        device,
        epochs : int,
        learning_rate : float,
        batch_size : int = None,
        silent_mode : bool = False,
        **kwargs
        ) -> NoReturn:

    model.to(device)
    optim = t.optim.Adam(model.parameters(),
                         lr = learning_rate)

    trackLoss = LossTracker()
    progressBar = SimpleProgressBar(epochs,
                                    silent_mode = silent_mode,
                                    length = 20)

    # use full dataset if no Batch Size specified
    if batch_size is None:
        batch_size = dataset.M
    else:
        batch_size = int(np.min((batch_size,dataset.M)))

    dataloader = DataLoader(dataset,
                            batch_size = batch_size,
                            shuffle = False,
                            )

    # Use try/except to catch SIGINT for early interuption
    try:
        for epoch in range(epochs):
            epoch_loss = 0.0
            for batch in dataloader:

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

            progressBar(epoch, epoch_loss)
            trackLoss(epoch_loss)

        # newline after complettion
        print('\n')

    except KeyboardInterrupt:
        print(f"\n\nPress Ctrl+C again to interrupt whole process")

    return trackLoss.history



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

    st_loss_history = fit(dataset = st_data,
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



    sc_loss_history = fit(dataset = sc_data,
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


