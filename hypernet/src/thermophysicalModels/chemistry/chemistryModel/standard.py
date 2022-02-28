import numpy as np

from hypernet.src.thermophysicalModels.chemistry.chemistryModel import Basic


class Standard(Basic):

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
        super(Standard, self).__init__(
            specieThermos,
            processFlags,
            reactionsList=reactionsList,
            *args,
            **kwargs
        )
        self.m = self.spTh[self.atom].specie.m

    # Methods
    ###########################################################################
    # Rates matrices ----------------------------------------------------------
    def K_(self, reac):
        labels = {
            'f': 'kf',
            'r': 'kr',
        }
        Ke = self.Ke_(self.processFlags['excit'], reac, labels) / self.m
        Kd = self.Kd_(self.processFlags['diss'], reac, labels) / self.m
        Kr = self.Kr_(self.processFlags['diss'], reac, labels) / self.m**2*2
        return Ke, Kd, Kr

    # Rates matrices derivatives ----------------------------------------------
    def dKdT_(self, reac):
        labels = {
            'f': 'dkfdT',
            'r': 'dkrdT',
        }
        dKedT = self.Ke_(self.processFlags['excit'], reac, labels) / self.m
        dKddT = self.Kd_(self.processFlags['diss'], reac, labels) / self.m
        dKrdT = self.Kr_(self.processFlags['diss'], reac, labels) / self.m**2*2
        return dKedT, dKddT, dKrdT

    # Porcesses matrices ------------------------------------------------------
    def Ke_(self, mask, reac, labels):
        '''Excit. & Relax. rates matrix'''

        # Construct Excit. & Relax. matrix
        K = np.zeros((self.nSpecies,self.nSpecies), dtype=np.float64)

        if mask:
            # Get excitation/relaxation rates
            reac = reac.loc[reac['reacIndex'].isin(self.processIndices['excit'])]

            # Fill matrix
            for i, row in reac.iterrows():
                l, r = row['indices']
                K[l,r] = K[l,r] + row[labels['f']]
                K[r,l] = K[r,l] + row[labels['r']]

            # Manipulate matrix
            K = -np.diag(np.sum(K, axis=1)) + np.transpose(K)

        return K

    def Kd_(self, mask, reac, labels):
        '''Dissociation rates matrix'''

        # Construct Dissociation matrix
        K = np.zeros((self.nSpecies,self.nSpecies), dtype=np.float64)

        if mask:
            # Get dissociation rates
            reac = reac.loc[reac['reacIndex'].isin(self.processIndices['diss'])]
            reac.sort_values(by=['indices'])

            # Extract forward rates
            rates = reac[labels['f']].to_numpy(dtype=np.float64)

            # Fill matrix
            K[:-1,:-1] = np.diag(-rates)
            K[ -1,:-1] = rates

        return K

    def Kr_(self, mask, reac, labels):
        '''Recombination rates matrix'''

        # Construct Recombination matrix
        K = np.zeros((self.nSpecies,1), dtype=np.float64)

        if mask:
            # Get Recombination rates
            reac = reac.loc[reac['reacIndex'].isin(self.processIndices['diss'])]
            reac.sort_values(by=['indices'])

            # Extract reverse rates
            rates = reac[labels['r']].to_numpy(dtype=np.float64)

            # Fill matrix
            K[:-1,0] = rates
            K[ -1,0] = -np.sum(rates)

        return K
