import os
import abc

from hypernet.src.general import utils
from hypernet.src.thermophysicalModels.chemistry.reactions.reactions import Reactions

import hypernet.database as db
kinetic_db = os.path.dirname(db.__file__) + '/air/kinetics/'


class Basic(object):

    # Initialization
    ###########################################################################
    def __init__(
        self,
        specieThermos,
        processFlags,
        reactionsList=None,
        *args,
        **kwargs
    ):
        # Thermodynamic specie properties
        self.spTh = specieThermos

        # Species
        self.nSpecies, self.specieIndices = self.species_params()
        for name, spTh in self.spTh.items():
            if spTh.specie.n_at > 1:
                self.molecule = name
            else:
                self.atom = name

        # Reactive processes
        self.processFlags = processFlags
        self.processIndices = {
            'diss':  [2],
            'excit': [5, 6]
        }

        # Reactions
        if reactionsList is None:
            rovib = self.spTh[self.molecule].specie.rovib
            self.reactionsList = kinetic_db + '/' + rovib['system'] + '_' \
                + rovib['PES'] + '/' + rovib['grouping'] + '/reactions.csv'
            utils.warning(
                "No specific `reactionsList` have been "
                "provided to the chemistry model."
            )
            utils.print_submain(">> The following one will be taken:")
            utils.print_submain(
                ">> `{}`".format(os.path.normpath(self.reactionsList))
            )
        else:
            self.reactionsList = reactionsList

        self.reactions = Reactions(
            self.spTh,
            self.reactionsList,
            self.processIndices,
            *args,
            **kwargs
        )

    # Properties
    ###########################################################################
    # Rates matrices ----------------------------------------------------------
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

    # Methods
    ###########################################################################
    # Update method -----------------------------------------------------------
    def update(self, T):
        reac = self.reactions.update(T)
        self.K = self.K_(reac)
        self.dKdT = self.dKdT_(reac)

    # Species parameters ------------------------------------------------------
    def species_params(self):
        nSpecies, specieIndices = 0, {}
        for name, spTh in self.spTh.items():
            if spTh.specie.n_at > 1:
                n = spTh.specie.n_bins
            else:
                n = 1
            specieIndices[name] = list(range(nSpecies, nSpecies+n))
            nSpecies += n
        return nSpecies, specieIndices

    # Rates matrices ----------------------------------------------------------
    @abc.abstractmethod
    def K_(self, reac):
        pass

    @abc.abstractmethod
    def dKdT_(self, reac):
        pass
