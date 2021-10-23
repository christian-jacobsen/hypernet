import numpy as np
import pandas as pd

from hypernet.src.general import const
from hypernet.src.general import utils

from hypernet.src.thermophysicalModels.specie import specie as specieModule
from hypernet.src.thermophysicalModels.specie import partitionFun as PFModule
from hypernet.src.thermophysicalModels.specie import thermo as thermoModule
from hypernet.src.thermophysicalModels.specie import equationOfState as EOSModule



class Standard(Basic):

    # Initialization
    ###########################################################################
    def __init__(
        self,
        specieThermos,
        reactionsList,
        processFlags,
        *args,
        **kwargs
    ):
        super(Standard, self).update(
            specieThermos,
            reactionRate,
            processIndeces,
            *args,
            **kwargs
        )

    # Master Equation matrices ------------------------------------------------
    def K_(self, reac):
        '''Get all the rates matrices for the ODE.'''
        m = self.spTh['O'].specie.m

        Ke = self.Ke_(self.processFlags['excit'], reac) / m
        Kd = self.Kd_(self.processFlags['diss'], reac) / m
        Kr = self.Kr_(self.processFlags['diss'], reac) / m**2 * 2
        return Ke, Kd, Kr

    def Ke_(self, mask, reac):
        '''Excit. & Relax. rates matrix'''

        # Construct Excit. & Relax. matrix
        K = np.zeros((self.nSpecies,self.nSpecies), dtype=np.float64)

        if mask:
            # Get excitation/relaxation rates
            reac = reac.loc[reac['reacIndex'].isin(self.reacIndeces['excit'])]

            # Fill matrix
            for i, row in reac.iterrows():
                l, r = row['indeces']
                K[l,r] += row['kf']
                K[r,l] += row['kr']

            # Manipulate matrix
            K = -np.diag(np.sum(K, axis=1)) + np.transpose(K)

        return K

    def Kd_(self, mask, reac):
        '''Dissociation rates matrix'''

        # Construct Dissociation matrix
        K = np.zeros((self.nSpecies,self.nSpecies), dtype=np.float64)

        if mask:
            # Get dissociation rates
            reac = reac.loc[reac['reacIndex'].isin(self.reacIndeces['diss'])]
            reac.sort_values(by=['indeces'])

            # Extract forward rates
            rates = reac['kf'].values

            # Fill matrix
            K[:-1,:-1] = np.diag(-rates)
            K[ -1,:-1] = rates

        return K

    def Kr_(self, mask, reac):
        '''Recombination rates matrix'''

        # Construct Recombination matrix
        K = np.zeros((self.nSpecies,1), dtype=np.float64)

        if mask:
            # Get Recombination rates
            reac = reac.loc[reac['reacIndex'].isin(self.reacIndeces['diss'])]
            reac.sort_values(by=['indeces'])

            # Extract reverse rates
            rates = reac['kr'].values

            # Fill matrix
            K[:-1,0] = rates
            K[ -1,0] = -np.sum(rates)

        return K

    # Master Equation matrices derivatives ------------------------------------
    def dKdT_(self, reac):
        '''Get all the rates matrices for the ODE.'''
        m = self.spTh['O'].specie.m

        dKedT = self.dKedT_(self.processFlags['excit'], reac) / m
        dKddT = self.dKddT_(self.processFlags['diss'], reac) / m
        dKrdT = self.dKrdT_(self.processFlags['diss'], reac) / m**2 * 2
        return dKedT, dKddT, dKrdT

    def dKedT_(self, mask, reac):
        '''Excit. & Relax. rates matrix'''

        # Construct Excit. & Relax. matrix
        dKdT = np.zeros((self.nSpecies,self.nSpecies), dtype=np.float64)

        if mask:
            # Get excitation/relaxation rates
            reac = reac.loc[reac['reacIndex'].isin(self.reacIndeces['excit'])]

            # Fill matrix
            for i, row in reac.iterrows():
                l, r = row['indeces']
                dKdT[l,r] += row['dkfdT']
                dKdT[r,l] += row['dkrdT']

            # Manipulate matrix
            dKdT = -np.diag(np.sum(dKdT, axis=1)) + np.transpose(dKdT)

        return dKdT

    def dKddT_(self, mask, reac):
        '''Dissociation rates matrix'''

        # Construct Dissociation matrix
        dKdT = np.zeros((self.nSpecies,self.nSpecies), dtype=np.float64)

        if mask:
            # Get dissociation rates
            reac = reac.loc[reac['reacIndex'].isin(self.reacIndeces['diss'])]
            reac.sort_values(by=['indeces'])

            # Extract reverse rates
            rates = reac['dkfdT'].values

            # Fill matrix
            dKdT[:-1,:-1] = np.diag(-rates)
            dKdT[ -1,:-1] = rates

        return dKdT

    def dKrdT_(self, mask, reac):
        '''Recombination rates matrix'''

        # Construct Recombination matrix
        dKdT = np.zeros((self.nSpecies,1), dtype=np.float64)

        if mask:
            # Obtain Recombination rates
            reac = reac.loc[reac['reacIndex'].isin(self.reacIndeces['diss'])]
            reac.sort_values(by=['indeces'])

            # Extract reverse rates
            rates = reac['dkrdT'].values

            # Fill matrix
            dKdT[:-1,0] = rates
            dKdT[ -1,0] = -np.sum(rates)

        return dKdT
