import abc

from hypernet.src.thermophysicalModels.chemistry.reactions import Reactions


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

    # Species details ---------------------------------------------------------
    def species(self):
        nSpecies, specieIndeces = 0, {}
        for name, spTh_ in self.spTh.values():
            if spTh_.specie.n_at > 1:
                n = spTh_.specie.n_bins
            else:
                n = 1
            specieIndeces[name] = list(range(nSpecies, nSpecies+n))
            nSpecies += n
        return nSpecies, specieIndeces

    # Rates matrices ----------------------------------------------------------
    @abc.abstractmethod
    def K_(self, reac):
        pass

    @abc.abstractmethod
    def dKdT_(self, reac):
        pass
