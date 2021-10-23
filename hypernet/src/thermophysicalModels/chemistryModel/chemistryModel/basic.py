import abc


import numpy as np
import pandas as pd

from hypernet.src.general import const
from hypernet.src.general import utils

from hypernet.src.thermophysicalModels.specie import specie as specieModule
from hypernet.src.thermophysicalModels.specie import partitionFun as PFModule
from hypernet.src.thermophysicalModels.specie import thermo as thermoModule
from hypernet.src.thermophysicalModels.specie import equationOfState as EOSModule




class Basic(object):

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
        # Thermodynamic specie properties
        self.spTh = specieThermos

        # Reactive processes
        self.processFlags = processFlags
        self.processIndeces = {
            'diss':  [2],
            'excit': [5, 6]
        }

        # Reactions
        self.reactionsList = reactionsList
        self.reactions = Reactions(
            self.spTh,
            self.reactionsList,
            self.processIndeces,
            *args,
            **kwargs
        )

        # Species
        self.nSpecies, self.specieIndeces = self.species()

    # Update method -----------------------------------------------------------
    def update(self, T):
        reac = self.reactions.update(T)
        self.K = self.K_(reac)
        self.dKdT = self.dKdT_(reac)

    # Properties --------------------------------------------------------------
    @property
    def K(self):
        return self._K
    @K.setter
    def K(self, value):
        self._K = value

    @property
    def dKdT(self):
        return self._dKdT
    @dKdT.setter
    def dKdT(self, value):
        self._dKdT = value

    # Reverse reaction rates --------------------------------------------------
    @abc.abstractmethod
    def K_(self, reac):
        pass

    @abc.abstractmethod
    def dKdT_(self, reac):
        pass
