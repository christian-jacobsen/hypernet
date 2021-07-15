import re
import os
import sys
import numpy as np
import pandas as pd
import inputs as inp
# Remove the following line if you install `hypernet` as a Python package -----
sys.path.append(inp.hypernet)

from tqdm import tqdm
from hypernet.src.general import utils


class DataGenerator(object):

    # Initialization
    ###########################################################################
    def __init__(
        self,
        T,
        species,
        rates_path,
        reac_files,
        reac_idx
    ):
        for name, specie in species.items():
            if specie.n_at > 1:
                self.n_bins = specie.n_bins

        # Data quantities =====================================================
        self.T_str = list(map(lambda x: str(int(x)), T))

        # Retrieving rates coefficients =======================================
        self.rates_path = rates_path
        self.react_file = reac_files

        # Building Master Equation matrices ===================================
        self.num_diss = int(self.n_bins)
        self.num_exc = int((self.n_bins-1)*self.n_bins/2)
        self.rates_coeff, self.rates_names = self.get_rates()
        self.reac_idx = np.concatenate(
            (np.full((self.num_diss,), reac_idx['diss']),
            np.full((self.num_exc,), reac_idx['exc']))
        )

    # Methods
    ###########################################################################
    # Training Data ===========================================================
    def training_data(self):

        x = np.zeros((len(self.T_str), 1), dtype=np.float64)
        y = np.zeros((len(self.T_str), self.num_exc + self.num_diss), dtype=np.float64)

        for i, T_i in enumerate(self.rates_coeff):
            x[i] = float(T_i)
            y[i,:] = np.hstack(tuple([K for name,K in self.rates_coeff[T_i].items()]))

        return x,y

    # Kinetic rates ===========================================================
    def get_rates(self):
        '''Retrieve rates coefficients from file.'''
        rates_coeff = dict()
        rates_name  = []
        utils.print_submain('Reading rates')
        for i, T_i in enumerate(tqdm(self.T_str)):
            if T_i in rates_coeff:
                continue
            rates_coeff[T_i] = dict()
            rates_name_i = []
            for reac in self.react_file:
                data_file = self.rates_path + 'T' + T_i + 'K/' + self.react_file[reac]
                rates_coeff[T_i][reac], reac_names = self.read_rates(data_file, reac)
                rates_name_i.extend(reac_names)
            if len(rates_name_i) > len(rates_name):
                rates_name = rates_name_i
            # Sum together 'exch' and 'inel' rates
            rates_coeff[T_i] = {
                'diss': rates_coeff[T_i]['diss'],
                'exc': rates_coeff[T_i]['exch']+rates_coeff[T_i]['inel'],
            }
        rates_name = rates_name[:self.num_exc + self.num_diss]
        return rates_coeff, rates_name

    def read_rates(self, file, reac_name):
        '''Read rates coefficients from file.'''
        shape = (1,self.num_diss) if reac_name == 'diss' else (1,self.num_exc)
        rates = np.zeros(shape, dtype=np.float64)

        infile = open(file, 'r')
        reac_eqn = []
        for i, line in enumerate(infile):
            _reac_eqn, rate_str = line.split(':')
            rate = float(rate_str.split(',')[0])
            bins_idx = re.findall(r"\((.*?)\)",_reac_eqn)
            if reac_name == 'diss':
                rates[0,int(bins_idx[0])-1] = rate
            else:
                rates[0,i] = rate
            reac_eqn.append(_reac_eqn)
        infile.close()

        return rates, reac_eqn
