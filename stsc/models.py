#!/usr/bin/env python3

import sys
import torch as t
import torch.nn as nn
from torch.nn.parameter import Parameter


import numpy as np
import pandas as pd

from typing import NoReturn, List, Tuple, Union, Collection
import logging

import os.path as osp


class ScModel(nn.Module):
    """ Model for singel cell data """

    def __init__(self,
                 n_genes : int,
                 n_celltypes : int,
                 device : t.device,
                 )->None:

        super().__init__()

        # Get dimensions from data
        self.K = n_celltypes
        self.G = n_genes

        # Define parameters to be estimated
        self.theta = Parameter(t.Tensor(self.G,self.K).to(device))
        self.R = t.Tensor(self.G,self.K).to(device)
        self.o = Parameter(t.Tensor(self.G,1).to(device))

        # Initialize parameters
        nn.init.normal_(self.o,
                        mean = 0.0,
                        std = 1.0)

        nn.init.normal_(self.theta,
                        mean = 0.0,
                        std = 1.0)



        # Functions to be used
        self.nb = t.distributions.NegativeBinomial
        self.softpl = nn.functional.softplus
        self.logsig = nn.functional.logsigmoid

    def _llnb(self,
              x : t.Tensor,
              meta : t.LongTensor,
              sf : t.Tensor,
             ) -> t.Tensor :

        """Log Likelihood for NB-model

        Returns the log likelihood for rates and logodds
        taken as a function of the observed counts.
        Assumes that single cell data is negative
        binomial distributed.

        Returns
        -------
        The log likelihood

        """


        log_unnormalized_prob = (sf*self.R[:,meta] * self.logsig(-self.o) +
                                 x * self.logsig(self.o))

        log_normalization = -t.lgamma(sf*self.R[:,meta] + x) + \
                             t.lgamma(1. + x) + \
                             t.lgamma(sf*self.R[:,meta])

        ll = t.sum(log_unnormalized_prob - log_normalization)

        return ll


    def forward(self,
                x : t.Tensor,
                meta : t.LongTensor,
                sf : t.Tensor,
                **kwargs,
                ) -> t.Tensor :
        """Forward pass during optimization"""

        # rates for each cell type
        self.R = self.softpl(self.theta)
        # get loss for current parameters
        self.loss = -self._llnb(x.transpose(1,0),
                                meta,
                                sf)

        return self.loss

    def __str__(self,):
        return f"sc_model"

class STModel(nn.Module):

    def __init__(self,
                 n_spots: int,
                 R : np.ndarray,
                 logits : np.ndarray,
                 device : t.device,
                 **kwargs,
                 )->None:

        super().__init__()

        # Get dimensions for from data
        self.S = n_spots
        self.G, self.K = R.shape
        self.Z = self.K + 1

        # Data from single cell estimates; Rates (R) and logits (o)
        self.R = t.tensor(R.astype(np.float32)).to(device)
        self.o = t.tensor(logits.astype(np.float32).reshape(-1,1)).to(device)

        # model specific parameters
        self.softpl = nn.functional.softplus
        self.lsig = nn.functional.logsigmoid
        self.sig = t.sigmoid

        # Learn noise from data
        self.eta = Parameter(t.tensor(np.zeros((self.G,1)).astype(np.float32)).to(device))
        nn.init.normal_(self.eta, mean = 0.0, std = 1.0)


        # un-normalized proportion in log space
        self.theta = Parameter(t.tensor(np.zeros((self.Z,self.S)).astype(np.float32)).to(device))
        nn.init.normal_(self.theta, mean = 0.0,std = 1.0)
        # gene bias in log space
        if not kwargs.get("freeze_beta",False):
            self.beta = Parameter(t.tensor(np.zeros((self.G,1)).astype(np.float32)).to(device))
            self.beta_trans = self.softpl
            nn.init.normal_(self.beta, mean = 0.0, std = 0.1)
        else:
            print("Using static beta_g")
            self.beta = t.tensor(np.ones((self.G,1)).astype(np.float32)).to(device)
            self.beta_trans = lambda x : x
        # un-normalized proportions
        self.v = t.tensor(np.zeros((self.Z,self.S)).astype(np.float32)).to(device)

        self.loss = t.tensor(0.0)
        self.model_ll = 0.0


    def noise_loss(self,
                  )-> t.Tensor:
        """Regularizing term for noise"""
        return -0.5*t.sum(t.pow(self.eta,2))

    def _llnb(self,
              x : t.Tensor,
              )->t.Tensor:
        """Log Likelihood function for standard model"""

        log_unnormalized_prob = self.r * self.lsig(-self.o) + \
                                x * self.lsig(self.o)

        log_normalization = -t.lgamma(self.r + x) + \
                             t.lgamma(1. + x) + \
                             t.lgamma(self.r)


        ll = t.sum(log_unnormalized_prob - log_normalization)

        self.ll = ll.item()

        return ll

    def _lfun(self,
              x : t.Tensor,
              )-> t.Tensor:
        """Loss Function

        Composed of the likelihood and prior of
        noise. Returns negative value of the above
        terms, to obtain a proper loss function.

        L(x) = -[LogLikelihood(x) + log(prior(noise))]

        Parameter
        --------
        x : t.tensor
            observed counts (n_genes x n_spots)

        """

        # log likelihood of observed count given model
        data_loss = self._llnb(x)
        # log of prior on noise elements
        noise_loss = self.noise_loss()

        return  - data_loss - noise_loss

    def __str__(self,
               )-> str:
        return f"st_model"

    def forward(self,
                x : t.tensor,
                gidx : t.tensor,
                **kwargs,
                ) -> t.tensor:

        """Forward pass"""

        self.gidx = gidx
        # proportion values
        self.v = self.softpl(self.theta)
        # noise values
        self.eps = self.softpl(self.eta)
        # account for gene specific bias and add noise
        self.Rhat = t.cat((t.mul(self.beta_trans(self.beta), self.R),self.eps),dim = 1)
        # combinde rates for all cell types
        self.r = t.einsum('gz,zs->gs',[self.Rhat,self.v[:,self.gidx]])
        # get loss for current parameters
        self.loss = self._lfun(x.transpose(1,0))

        return self.loss



