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

        # Reactive processes
        self.processFlags = processFlags
        self.processIndices = {
            'diss':  [2],
            'excit': [5, 6]
        }

        # Reactions
        self.reactionsList = reactionsList
        if self.reactionsList is None:
            for spTh_ in self.spTh.values():
                if spTh_.specie.n_at > 1:
                    molecule = spTh_.specie
                break
            self.reactionsList = kinetic_db + '/' \
                + molecule.rovib['system'] + '_' + molecule.rovib['PES'] \
                + '/' + molecule.rovib['grouping'] + '/reactions.csv'
            utils.warning("No specific `reactionsList` have been provided to "
                "the chemistry model. The following one will be taken:\n"
                "`{}`.".format(os.path.normpath(self.reactionsList)))

        self.reactions = Reactions(
            self.spTh,
            self.reactionsList,
            self.processIndices,
            *args,
            **kwargs
        )

        # Species
        self.nSpecies, self.specieIndices = self.species()

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
        print(reac)
        self.K = self.K_(reac)
        self.dKdT = self.dKdT_(reac)

    # Species details ---------------------------------------------------------
    def species(self):
        nSpecies, specieIndices = 0, {}
        for name, spTh_ in self.spTh.items():
            if spTh_.specie.n_at > 1:
                n = spTh_.specie.n_bins
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
