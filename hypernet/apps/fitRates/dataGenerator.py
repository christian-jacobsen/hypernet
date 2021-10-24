import os
import re
import numpy as np
import pandas as pd

from hypernet.src.general import utils

import hypernet.database as db
kinetic_db = os.path.dirname(db.__file__) + '/air/kinetics/'


class DataGenerator(object):

    # Initialization
    ###########################################################################
    def __init__(
        self,
        species,
        T,
        reacReader,
        reacWriter,
        verbose=True
    ):
        # Verbosity
        self.verbose = verbose

        # Species details
        self.species = species
        for name, sp in species.items():
            if sp.n_at > 1:
                self.n_bins = sp.n_bins

        # Temperatures
        self.TStr = [ str(int(i)) for i in T ]

        # Reaction reader/writer
        self.reacReader = reacReader
        self.reacWriter = reacWriter

        # Define number of reactions
        self.n_diss = int(self.n_bins)
        self.n_excit = int((self.n_bins-1)*self.n_bins/2)
        self.n_reac = self.n_excit*2 + self.n_diss

        # Reactions database header
        self.header = ['description', 'reactStr', 'reacIndex']

        # Reactions database
        self.rates, self.reacDB = self.get_rates()

    # Methods
    ###########################################################################
    # Training Data -----------------------------------------------------------
    def training_data(self):

        x = np.zeros((len(self.TStr), 1), dtype=np.float64)
        y = np.zeros((len(self.TStr), self.n_reac), dtype=np.float64)

        for i, T_i in enumerate(self.rates):
            x[i,0] = float(T_i)
            y[i,:] = np.hstack(
                tuple([K for name,K in self.rates[T_i].items()])
            )

        return x, y

    # Kinetic rates -----------------------------------------------------------
    def get_rates(self):
        '''Retrieve rates coefficients from file.'''
        rates = dict()
        for i, T_i in enumerate(self.TStr):
            if T_i in rates:
                continue
            rates[T_i] = dict()
            utils.print_main('> Temp.: %5d K' % int(T_i), verbose=self.verbose)
            reacDB = pd.DataFrame(columns=self.header)
            for process, file in self.reacReader['files'].items():
                path = kinetic_db+self.reacReader['path']+'T'+T_i+'K/'+file
                rates[T_i][process], reacDB_ = self.read_rates(
                    path, process
                )
                reacDB = pd.concat(
                    [reacDB, pd.DataFrame.from_dict(reacDB_)], ignore_index=True
                )
        return rates, reacDB

    def read_rates(self, path, process):
        '''Read rates coefficients from file.'''
        utils.print_submain(
            'Reading `{}` rates from file'.format(process),
            verbose=self.verbose
        )

        # Initialize dataframe
        reacDB = {k: [] for k in self.header}

        # Initialize rates vector
        shape = (1, self.n_diss) if process == 'diss' else (1, self.n_excit)
        rates = np.zeros(shape, dtype=np.float64)

        # Open/read file
        infile = open(path, 'r')
        for i, line in enumerate(infile):
            # Extract reaction rate and index
            if self.reacReader['file_version'] == 'old':
                eqn, rate = line.split(':')
                rate = float(rate.split(',')[0]) * 1.e-6
                bins_idx = re.findall(r"\((.*?)\)", eqn)
                rate_idx = int(bins_idx[0])-1 if process == 'diss' else i
            else:
                if i == 0:
                    continue
                *bins_idx, rate = line.split(',')
                rate = float(rate) * 1.e-6
                rate_idx = int(bins_idx[0])-1 if process == 'diss' else i-1
            # Fill rates vector
            rates[0,rate_idx] = rate
            # Update dataframe
            reacDB = self.updateDB(reacDB, bins_idx, process)

        infile.close()
        return rates, reacDB

    def updateDB(self, reacDB, bins_idx, process):
        reacDB['description'].append(
            self.reacWriter['type']+self.reacWriter['rate']
        )
        reacDB['reactStr'].append(
            self.compose_reacStr(bins_idx, process)
        )
        reacDB['reacIndex'].append(
            self.reacWriter['indices'][process]
        )
        return reacDB

    def compose_reacStr(self, bins_idx, process):
        n_at = 0
        for name, sp in self.species.items():
            n_at += sp.n_at
            if sp.n_at > 1:
                molecule = name
            else:
                atom = name

        lhs = molecule + '(' + bins_idx[0] + ')+' + atom
        if process == 'diss':
            rhs = str(n_at) + atom
        else:
            rhs = molecule + '(' + bins_idx[1] + ')+' + atom

        return lhs + '=' + rhs
